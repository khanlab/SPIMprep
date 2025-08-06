import json
import os
import re
from glob import glob
from itertools import product
from pathlib import Path

import tifffile
import xmltodict
from snakemake.io import glob_wildcards

tif_files = glob(f"{snakemake.input.ome_dir}/*.tif")

# check first tif file to see if it is zstack or not:
if "xyz-Table Z" in Path(tif_files[0]).name:
    is_zstack = False
else:
    is_zstack = True

# now check to see if there is just single tile
if " x " in Path(tif_files[0]).name:
    is_tiled = True
else:
    is_tiled = False

if is_tiled:
    if is_zstack:
        in_tif_pattern = os.path.join(
            snakemake.input.ome_dir,
            snakemake.config["import_blaze"]["raw_tif_pattern_zstack"],
        )
    else:
        in_tif_pattern = os.path.join(
            snakemake.input.ome_dir, snakemake.config["import_blaze"]["raw_tif_pattern"]
        )
else:
    in_tif_pattern = os.path.join(
        snakemake.input.ome_dir,
        snakemake.config["import_blaze"]["raw_tif_pattern_notile"],
    )

# add a wildcard constraint to ensure no
# subfolders get parsed (ie don't match anything with / in it):
prefix_constraint = r"[^/]+"
in_tif_pattern_constrained = in_tif_pattern.replace(
    "{prefix}", f"{{prefix,{prefix_constraint}}}"
)


# parse the filenames to get number of channels, tiles etc..
if is_tiled:
    if is_zstack:
        prefix, tilex, tiley, channel = glob_wildcards(
            in_tif_pattern_constrained, tif_files
        )
    else:
        prefix, tilex, tiley, channel, zslice = glob_wildcards(
            in_tif_pattern_constrained, tif_files
        )
else:
    prefix, channel, zslice = glob_wildcards(in_tif_pattern_constrained, tif_files)

if is_tiled:
    tiles_x = sorted(list(set(tilex)))
    tiles_y = sorted(list(set(tiley)))

channels = sorted(list(set(channel)))
prefixes = sorted(list(set(prefix)))


if not is_zstack:
    zslices = sorted(list(set(zslice)))

if is_tiled:
    if is_zstack:
        in_tif = in_tif_pattern.format(
            tilex=tiles_x[0], tiley=tiles_y[0], prefix=prefixes[0], channel=channels[0]
        )
    else:
        # read in series metadata from first file
        in_tif = in_tif_pattern.format(
            tilex=tiles_x[0],
            tiley=tiles_y[0],
            prefix=prefixes[0],
            channel=channels[0],
            zslice=zslices[0],
        )
else:
    in_tif = in_tif_pattern.format(
        prefix=prefixes[0], channel=channels[0], zslice=zslices[0]
    )


raw_tif = tifffile.TiffFile(in_tif, mode="r")

axes = raw_tif.series[0].get_axes()
shape = raw_tif.series[0].get_shape()


ome_dict = xmltodict.parse(raw_tif.ome_metadata)
physical_size_x = ome_dict["OME"]["Image"]["Pixels"]["@PhysicalSizeX"]
physical_size_y = ome_dict["OME"]["Image"]["Pixels"]["@PhysicalSizeY"]
physical_size_z = ome_dict["OME"]["Image"]["Pixels"]["@PhysicalSizeZ"]
custom_metadata = ome_dict["OME"]["Image"]["ca:CustomAttributes"]


# read tile configuration from the microscope metadata
if axes == "CZYX":
    if is_tiled:
        if is_zstack:
            tile_config_pattern = r"Blaze\[(?P<tilex>[0-9]+) x (?P<tiley>[0-9]+)\]_C(?P<channel>[0-9]+).ome.tif;;\((?P<x>\S+), (?P<y>\S+),(?P<chan>\S+)\)"
        else:
            tile_config_pattern = r"Blaze\[(?P<tilex>[0-9]+) x (?P<tiley>[0-9]+)\]_C(?P<channel>[0-9]+)_xyz-Table Z(?P<zslice>[0-9]+).ome.tif;;\((?P<x>\S+), (?P<y>\S+),(?P<chan>\S+), (?P<z>\S+)\)"
    else:
        print("single tile, axes=CZYX")

elif axes == "ZYX":
    if is_tiled:
        if is_zstack:
            tile_config_pattern = r"Blaze\[(?P<tilex>[0-9]+) x (?P<tiley>[0-9]+)\]_C(?P<channel>[0-9]+).ome.tif;;\((?P<x>\S+), (?P<y>\S+)\)"
        else:
            tile_config_pattern = r"Blaze\[(?P<tilex>[0-9]+) x (?P<tiley>[0-9]+)\]_C(?P<channel>[0-9]+)_xyz-Table Z(?P<zslice>[0-9]+).ome.tif;;\((?P<x>\S+), (?P<y>\S+), (?P<z>\S+)\)"
    else:
        print("single tile, axes=ZYX")

if is_tiled:
    tile_pattern = re.compile(tile_config_pattern)

    # put it in 3 maps, one for each coord, indexed by tilex, tiley, channel, and aslice
    map_x = dict()
    map_y = dict()
    map_z = dict()

    map_tiles_to_chunk = dict()
    chunks = []
    for chunk, (tilex, tiley) in enumerate(product(tiles_x, tiles_y)):
        map_tiles_to_chunk[tilex + tiley] = chunk
        chunks.append(chunk)

    for line in custom_metadata["TileConfiguration"]["@TileConfiguration"].split("  ")[
        1:
    ]:
        d = re.search(tile_pattern, line).groupdict()
        chunk = map_tiles_to_chunk[
            d["tilex"] + d["tiley"]
        ]  # want the key to have chunk instad of tilex,tiley, so map to that first

        if is_zstack:
            key = f"tile-{chunk}_chan-{d['channel']}_z-0000"
        else:
            # key is:  tile-{chunk}_chan-{channel}_z-{zslice}
            key = f"tile-{chunk}_chan-{d['channel']}_z-{d['zslice']}"

        map_x[key] = float(d["x"])
        map_y[key] = float(d["y"])
        if is_zstack:
            map_z[key] = float(0)
        else:
            map_z[key] = float(d["z"])


metadata = {}
if is_tiled:
    metadata["tiles_x"] = tiles_x
    metadata["tiles_y"] = tiles_y
    metadata["chunks"] = chunks
    metadata["lookup_tile_offset_x"] = map_x
    metadata["lookup_tile_offset_y"] = map_y
    metadata["lookup_tile_offset_z"] = map_z

metadata["channels"] = channels

if not is_zstack:
    metadata["zslices"] = zslices
metadata["prefixes"] = prefixes
metadata["axes"] = axes
metadata["shape"] = shape
metadata["physical_size_x"] = float(physical_size_x)
metadata["physical_size_y"] = float(physical_size_y)
metadata["physical_size_z"] = float(physical_size_z)
metadata["ome_full_metadata"] = ome_dict
metadata["PixelSize"] = [
    metadata["physical_size_z"] / 1000.0,
    metadata["physical_size_y"] / 1000.0,
    metadata["physical_size_x"] / 1000.0,
]  # zyx since OME-Zarr is ZYX
metadata["PixelSizeUnits"] = "mm"
metadata["is_zstack"] = is_zstack
metadata["is_tiled"] = is_tiled

# write metadata to json
with open(snakemake.output.metadata_json, "w") as fp:
    json.dump(metadata, fp, indent=4)
