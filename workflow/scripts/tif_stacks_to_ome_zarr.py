import json

import dask.array as da
from dask.array.image import imread as dask_imread
from dask.diagnostics import ProgressBar
from zarrnii import ZarrNii

in_tif_glob = snakemake.params.in_tif_glob
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
for i, stain in enumerate(stains):

    tif_glob = in_tif_glob.format(stain=stain)
    darr_list.append(dask_imread(tif_glob).rechunk(rechunk_size))

    # append to omero metadata
    channel_metadata = {
        key: val
        for key, val in snakemake.config["ome_zarr"]["omero_metadata"]["channels"][
            "defaults"
        ].items()
    }
    channel_name = stains[i]
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

# orientation="SAR" is ZYX order: z↑=Superior, y↑=Anterior, x↑=Right
# (equivalent to xyz_orientation="RAS" in new zarrnii >=0.2.0)
# NOTE: x/y directions assumed from LifeCanvas convention — verify
# visually in a NIfTI viewer after first run. z-axis (Superior) is
# correct for data where ascending TIF index = moving superiorly.
znimg = ZarrNii.from_darr(
    darr_channels,
    orientation="SAR",
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
    znimg.to_ome_zarr(out_zarr, version="0.4")