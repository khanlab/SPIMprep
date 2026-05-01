import os
import re
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import zarr
from tifffile import TiffFile

# ----------------------------
# User settings
# ----------------------------
# INPUT_DIR = Path("/nfs/trident3/lightsheet/prado/mouse_app_lecanemab_ki3_batch3/raw/tif_4x/A_AS38F4")
# OUTPUT_ZARR = Path("output_dataset.zarr")

INPUT_DIR = Path(snakemake.input.tif_dir)
OUTPUT_ZARR = Path(snakemake.output.zarr)

MAX_WORKERS = 16
Z_CHUNK = 16
DATASET_NAME = "0"

TEST_N_ROWS = 1000
TEST_N_COLS = 1000

# TODO: get this from config..
FNAME_RE = re.compile(r"Blaze\[(\d+)\s*x\s*(\d+)\]_C(\d+)\.ome\.tif$")


def parse_file_info(path: Path):
    m = FNAME_RE.search(path.name)
    if not m:
        return None
    row = int(m.group(1))
    col = int(m.group(2))
    ch = int(m.group(3))
    return row, col, ch


def validate_and_get_shape(path: Path):
    with TiffFile(path) as tif:
        pages = tif.pages
        if len(pages) == 0:
            raise ValueError(f"No TIFF pages in {path}")

        p0 = pages[0]

        if p0.compression != 1:
            raise ValueError(
                f"{path}: expected uncompressed TIFF, got compression={p0.compression}"
            )
        if p0.is_tiled:
            raise ValueError(f"{path}: expected striped TIFF, got tiled TIFF")

        first = p0.dataoffsets[0]
        last_end = p0.dataoffsets[-1] + p0.databytecounts[-1]
        expected = sum(p0.databytecounts)
        contiguous = (last_end - first) == expected
        if not contiguous:
            raise ValueError(f"{path}: page payload is not contiguous")

        z = len(pages)
        y, x = p0.shape
        dtype = np.dtype(p0.dtype)
        page_nbytes = expected

    return z, y, x, dtype, page_nbytes


def worker_write_one(
    tif_path_str: str,
    output_zarr_str: str,
    dataset_name: str,
    row: int,
    col: int,
    n_cols: int,
    ch: int,
    z_chunk: int,
):
    tif_path = Path(tif_path_str)
    tile = row * n_cols + col

    with TiffFile(tif_path) as tif:
        pages = tif.pages
        z = len(pages)
        p0 = pages[0]
        y, x = p0.shape
        dtype = np.dtype(p0.dtype)
        page_nbytes = sum(p0.databytecounts)

        if p0.compression != 1:
            raise ValueError(f"{tif_path}: expected uncompressed TIFF")
        if p0.is_tiled:
            raise ValueError(f"{tif_path}: expected non-tiled TIFF")

        first = p0.dataoffsets[0]
        last_end = p0.dataoffsets[-1] + p0.databytecounts[-1]
        if (last_end - first) != page_nbytes:
            raise ValueError(f"{tif_path}: page payload not contiguous")

        root = zarr.open_group(str(output_zarr_str), mode="r+")
        arr = root[dataset_name]

        if tile >= arr.shape[0] or ch >= arr.shape[1]:
            raise IndexError(
                f"Target index out of bounds for {tif_path.name}: "
                f"(tile={tile}, ch={ch}) vs array shape {arr.shape}"
            )

        t0_total = time.time()

        with open(tif_path, "rb", buffering=16 * 1024 * 1024) as fh:
            block = np.empty((z_chunk, y, x), dtype=dtype)

            for z0 in range(0, z, z_chunk):
                z1 = min(z0 + z_chunk, z)
                nz = z1 - z0

                for i, zi in enumerate(range(z0, z1)):
                    p = pages[zi]
                    fh.seek(p.dataoffsets[0])
                    raw = fh.read(page_nbytes)
                    block[i, :, :] = np.frombuffer(raw, dtype=dtype).reshape(y, x)

                arr[tile, ch, z0:z1, :, :] = block[:nz]

        dt = time.time() - t0_total
        return {
            "file": tif_path.name,
            "row": row,
            "col": col,
            "tile": tile,
            "ch": ch,
            "seconds": dt,
        }


def main():
    tif_files = sorted(INPUT_DIR.glob("*.ome.tif"))
    if not tif_files:
        raise SystemExit(f"No .ome.tif files found in {INPUT_DIR}")

    parsed = []
    for path in tif_files:
        info = parse_file_info(path)
        if info is None:
            print(f"Skipping unrecognized filename: {path.name}")
            continue
        row, col, ch = info
        parsed.append((path, row, col, ch))

    parsed = [
        (path, row, col, ch)
        for path, row, col, ch in parsed
        if row < TEST_N_ROWS and col < TEST_N_COLS
    ]

    if not parsed:
        raise SystemExit("No TIFFs matched the Blaze[row x col]_Cnn.ome.tif pattern")

    max_row = max(r for _, r, _, _ in parsed)
    max_col = max(c for _, _, c, _ in parsed)
    max_ch = max(ch for _, _, _, ch in parsed)

    n_rows = max_row + 1
    n_cols = max_col + 1
    n_channels = max_ch + 1
    n_tiles = n_rows * n_cols

    z, y, x, dtype, page_nbytes = validate_and_get_shape(parsed[0][0])

    print(f"Found {len(parsed)} TIFFs")
    print(f"Grid: rows={n_rows}, cols={n_cols}, tiles={n_tiles}, channels={n_channels}")
    print(f"Per-file stack: Z={z}, Y={y}, X={x}, dtype={dtype}")
    print(f"Using MAX_WORKERS={MAX_WORKERS}, Z_CHUNK={Z_CHUNK}")

    root = zarr.open_group(str(OUTPUT_ZARR), mode="w")
    arr = root.create_array(
        DATASET_NAME,
        shape=(n_tiles, n_channels, z, y, x),
        chunks=(1, 1, Z_CHUNK, y, x),
        dtype=dtype,
        compressors=None,
        overwrite=True,
    )

    root.attrs["input_dir"] = str(INPUT_DIR)
    root.attrs["tile_layout"] = {
        "n_rows": n_rows,
        "n_cols": n_cols,
        "tile_index_formula": "tile = row * n_cols + col",
    }

    t0 = time.time()
    failures = []

    with ProcessPoolExecutor(max_workers=MAX_WORKERS) as ex:
        futures = [
            ex.submit(
                worker_write_one,
                str(path),
                str(OUTPUT_ZARR),
                DATASET_NAME,
                row,
                col,
                n_cols,
                ch,
                Z_CHUNK,
            )
            for path, row, col, ch in parsed
        ]

        for fut in as_completed(futures):
            try:
                result = fut.result()
                print(
                    f"DONE {result['file']} "
                    f"(tile={result['tile']}, r={result['row']}, c={result['col']}, ch={result['ch']}) "
                    f"in {result['seconds']:.2f}s"
                )
            except Exception as e:
                failures.append(repr(e))
                print(f"FAILED: {type(e).__name__}: {e!r}")

    total_dt = time.time() - t0
    print(f"TOTAL elapsed: {total_dt:.2f}s")

    if failures:
        raise RuntimeError(f"{len(failures)} file(s) failed:\n" + "\n".join(failures))


if __name__ == "__main__":
    main()
