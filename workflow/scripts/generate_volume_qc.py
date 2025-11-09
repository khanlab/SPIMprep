import json
import math
import shutil
from distutils.dir_util import copy_tree
from pathlib import Path

import dask.array as da
import zarr
from lib.cloud_io import get_fsspec, is_remote
from upath import UPath as Path

from zarrnii import ZarrNii

# directory containing the volume rendering files
resource_dir = snakemake.output.resources
# where html file should be written
html_dest = snakemake.output.html


# move volume renderer into the subjects directory
copy_tree(snakemake.input.vol_viewer_dir, resource_dir)
shutil.move(Path(resource_dir) / "volRender.html", html_dest)

uri = snakemake.params.uri
if is_remote(uri):
    fs_args = {
        "storage_provider_settings": snakemake.params.storage_provider_settings,
        "creds": snakemake.input.creds,
    }
else:
    fs_args = {}

fs = get_fsspec(uri, **fs_args)

if Path(uri).suffix == ".zip":
    store = zarr.storage.ZipStore(Path(uri).path, dimension_separator="/", mode="r")
else:
    store = zarr.storage.FSStore(
        Path(uri).path, fs=fs, dimension_separator="/", mode="r"
    )
darr = da.from_zarr(store, component="/5")

# Get most downsampled ome-zarr image
ds_z = ZarrNii.from_darr(darr)
z_length = ds_z.darr.shape[1]

# downsample it so it has at most 100 slices and ast least 50 slices in z-direction
if z_length > 50:
    downsample_factor = math.floor(z_length / 50)
else:
    downsample_factor = 1
ds_z = ds_z.downsample(along_z=downsample_factor)

# Write it to a JSON for js script to read
with open(Path(resource_dir) / "volumeData.json", "w") as f:
    json_data = json.dumps(ds_z.darr.compute().tolist())
    f.write(json_data)
