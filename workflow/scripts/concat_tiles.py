from pathlib import Path

import dask.array as da
from dask.diagnostics import ProgressBar
from ngff_zarr.from_ngff_zarr import from_ngff_zarr

# Directory containing all sub-Zarr datasets
tiles_dir = Path(snakemake.input.tiles_dir)

arrays = []

for subdir in sorted(tiles_dir.iterdir()):
    zarr_path = subdir / "0"
    if zarr_path.is_dir():
        # Read the multiscale dataset (assumes first scale is highest resolution)
        multiscale = from_ngff_zarr(
            str(zarr_path)
        )  # Ensure string path if required by the library
        highest_res = multiscale.images[0].data  # This is a Dask array
        print(highest_res)
        arrays.append(highest_res)

# Concatenate all volumes along axis 0 (e.g., samples/timepoints/tiles/etc.)
stacked = da.concatenate(arrays, axis=0)

print(stacked)
# Save to new Zarr store (no metadata needed)

with ProgressBar():
    stacked.to_zarr(
        snakemake.output.zarr, overwrite=True
    )  # Convert to string if API expects string
