import json
import re
from itertools import product

import tifffile
import xmltodict
from snakemake.io import glob_wildcards

metadata = {}
metadata["physical_size_x"] = snakemake.params.physical_size_x_um
metadata["physical_size_y"] = snakemake.params.physical_size_y_um
metadata["physical_size_z"] = snakemake.params.physical_size_z_um
metadata["PixelSize"] = [
    float(metadata["physical_size_z"] / 1000.0),
    float(metadata["physical_size_y"] / 1000.0),
    float(metadata["physical_size_x"] / 1000.0),
]  # zyx since OME-Zarr is ZYX
metadata["PixelSizeUnits"] = "mm"

# write metadata to json
with open(snakemake.output.metadata_json, "w") as fp:
    json.dump(metadata, fp, indent=4)
