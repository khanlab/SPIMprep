import json

import dask.array as da
import zarr
from dask.diagnostics import ProgressBar
from lib.cloud_io import get_fsspec, is_remote
from upath import UPath as Path
from zarrnii import ZarrNii

in_zarr = snakemake.input.zarr


metadata_json = snakemake.input.metadata_json
downsampling = snakemake.params.downsampling
max_layer = (
    snakemake.params.max_downsampling_layers
)  # number of downsamplings by 2 to include in zarr
rechunk_size = snakemake.params.rechunk_size
out_zarr = snakemake.output.zarr
stains = snakemake.params.stains
scaling_method = snakemake.params.scaling_method

uri = snakemake.params.uri

# prepare metadata for ome-zarr
with open(metadata_json) as fp:
    metadata = json.load(fp)


# init omero metadata
omero = {
    key: val
    for key, val in snakemake.config["ome_zarr"]["omero_metadata"]["defaults"].items()
}
omero["channels"] = []


darr_list = []
for zarr_i in range(len(snakemake.input.zarr)):
    # open zarr to get group name
    in_zarr = snakemake.input.zarr[zarr_i]
    zi = zarr.open(in_zarr)
    group_name = [g for g in zi.group_keys()][0]

    darr_list.append(
        da.from_zarr(in_zarr, component=f"{group_name}/s0", chunks=rechunk_size)
    )

    # append to omero metadata
    channel_metadata = {
        key: val
        for key, val in snakemake.config["ome_zarr"]["omero_metadata"]["channels"][
            "defaults"
        ].items()
    }
    channel_name = stains[zarr_i]
    channel_metadata["label"] = channel_name
    default_color = snakemake.config["ome_zarr"]["omero_metadata"]["channels"][
        "default_color"
    ]
    color = snakemake.config["ome_zarr"]["omero_metadata"]["channels"][
        "color_mapping"
    ].get(channel_name, default_color)
    channel_metadata["color"] = color
    omero["channels"].append(channel_metadata)


darr_channels = da.stack(darr_list)


znimg = ZarrNii.from_darr(
    darr_channels,
    orientation="IPR",
    xyz_orientation="RPI",
    axes_order="ZYX",
    spacing=(
        float(metadata["physical_size_z"]) * float(downsampling) / 1000.0,
        float(metadata["physical_size_y"]) * float(downsampling) / 1000.0,
        float(metadata["physical_size_x"]) * float(downsampling) / 1000.0,
    ),
    unit="millimeter",
    omero=omero,
)

with ProgressBar():
    znimg.to_ome_zarr(snakemake.output[0], version="0.4")
