import logging
import sys
from pathlib import Path

import h5py
import hdf5plugin
import zarr

# ----------------------------
# Logging setup
# ----------------------------


use_stdout = True  # toggle this

if use_stdout:
    logging.basicConfig(
        stream=sys.stdout,
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
    )
else:
    log_path = snakemake.log[0] if isinstance(snakemake.log, list) else snakemake.log
    logging.basicConfig(
        filename=log_path,
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
    )

logger = logging.getLogger(__name__)

logger.info("Starting Imaris (HDF5) → Zarr conversion")


# ----------------------------
# Inputs
# ----------------------------
src_path = snakemake.input.ims
dst_path = snakemake.output.zarr
chan = snakemake.params.channel

h5_key = f"DataSet/ResolutionLevel 0/TimePoint 0/Channel {chan}/Data"

logger.info(f"Input file: {src_path}")
logger.info(f"Output path: {dst_path}")
logger.info(f"Channel: {chan}")
logger.info(f"HDF5 key: {h5_key}")


# ----------------------------
# Open source
# ----------------------------
source = h5py.File(src_path, mode="r")
src = source[h5_key]

logger.info(f"Source shape: {src.shape}")
logger.info(f"Source dtype: {src.dtype}")
logger.info(f"Source chunks: {src.chunks}")


# ----------------------------
# Store selection
# ----------------------------
if Path(dst_path).suffix == ".zip":
    logger.info("Using ZipStore")
    store = zarr.storage.ZipStore(dst_path, mode="w")
else:
    logger.info("Using LocalStore")
    store = zarr.storage.LocalStore(dst_path)

root = zarr.create_group(store=store, zarr_format=3)


# ----------------------------
# Chunking
# ----------------------------
chunks = src.chunks if src.chunks is not None else (256, 256, 256)

logger.info(f"Destination chunks: {chunks}")


# ----------------------------
# Create array
# ----------------------------
z = root.create_array(
    name="data",
    shape=src.shape,
    dtype=src.dtype,
    chunks=chunks,
    compressors=(),  # no compression
)

logger.info("Created destination Zarr array")


# ----------------------------
# Chunked copy
# ----------------------------
cz, cy, cx = z.chunks
sz, sy, sx = src.shape

total_blocks = ((sz + cz - 1) // cz) * ((sy + cy - 1) // cy) * ((sx + cx - 1) // cx)

logger.info(f"Total blocks to copy: {total_blocks}")

block_idx = 0

for z0 in range(0, sz, cz):
    z1 = min(z0 + cz, sz)
    for y0 in range(0, sy, cy):
        y1 = min(y0 + cy, sy)
        for x0 in range(0, sx, cx):
            x1 = min(x0 + cx, sx)

            block = src[z0:z1, y0:y1, x0:x1]
            z[z0:z1, y0:y1, x0:x1] = block

            block_idx += 1
            if block_idx % 1000 == 0 or block_idx == total_blocks:
                logger.info(f"Copied block {block_idx}/{total_blocks}")


# ----------------------------
# Cleanup
# ----------------------------
source.close()
logger.info("Closed source file")

# Ensure ZipStore flushes properly
if hasattr(store, "close"):
    store.close()
    logger.info("Closed Zarr store")

logger.info("Conversion complete")
