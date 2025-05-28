import tifffile
import os
import json
import dask.array as da
from itertools import product
from dask.diagnostics import ProgressBar
import dask
from glob import glob

def replace_square_brackets(pattern):
    """replace all [ and ] in the string (have to use 
    intermediate variable to avoid conflicts)"""
    pattern = pattern.replace('[','##LEFTBRACKET##')
    pattern = pattern.replace(']','##RIGHTBRACKET##')
    pattern = pattern.replace('##LEFTBRACKET##','[[]')
    pattern = pattern.replace('##RIGHTBRACKET##','[]]')
    return pattern


#read metadata json
with open(snakemake.input.metadata_json) as fp:
    metadata = json.load(fp)


is_zstack = metadata['is_zstack']
is_tiled = metadata['is_tiled']

if is_tiled:
    if is_zstack:
        in_tif_pattern = os.path.join(snakemake.input.ome_dir,snakemake.config["import_blaze"]["raw_tif_pattern_zstack"])
    else:
        in_tif_pattern = os.path.join(snakemake.input.ome_dir,snakemake.config["import_blaze"]["raw_tif_pattern"])
else:
    in_tif_pattern = os.path.join(snakemake.input.ome_dir,snakemake.config["import_blaze"]["raw_tif_pattern_notile"])



#use tif pattern but replace the [ and ] with [[] and []] so glob doesn't choke
in_tif_glob = replace_square_brackets(str(in_tif_pattern))



if is_tiled:

    size_tiles=len(metadata['tiles_x'])*len(metadata['tiles_y'])

    #now get the first channel and first zslice tif 
    tiles=[]
    for i_tile,(tilex,tiley) in enumerate(product(metadata['tiles_x'],metadata['tiles_y'])):
        print(f'tile {tilex}x{tiley}, {i_tile}') 
        channel=metadata['channels'][0]

        if is_zstack:
            tif_file = in_tif_pattern.format(tilex=tilex,tiley=tiley,prefix=metadata['prefixes'][0],channel=channel)
            
            in_store = tifffile.imread(tif_file, aszarr=True)  # this is a group

            tiles.append(da.from_zarr(in_store))

            print(tiles[-1].shape)
        else:
            tif_file = sorted(glob(in_tif_glob.format(tilex=tilex,tiley=tiley,prefix=metadata['prefixes'][0],channel=channel,zslice='*')))[0]
            in_store = tifffile.imread(tif_file, aszarr=True)  # this is a group

            tiles.append(da.from_zarr(in_store))


    print(tiles)


    #now we have list of tiles, each a dask array
    #stack them up to get our final array
    darr = da.stack(tiles)

else:
    #single tile, zslices:

    tif_file = sorted(glob(in_tif_glob.format(tilex=tilex,tiley=tiley,prefix=metadata['prefixes'][0],channel=channel,zslice='*')))[0]
    in_store = tifffile.imread(tif_file, aszarr=True)  # this is a group

    darr = da.from_zarr(in_store)




print(darr.shape)
print(darr.chunksize)

#rescale intensities, and recast 
darr = darr * snakemake.params.intensity_rescaling
darr = darr.astype('uint16')
darr = darr.rechunk((1,1,1,darr.shape[-2],darr.shape[-1]))

print(darr.chunksize)
#now we can do the computation itself, storing to zarr
print('writing images to zarr with dask')
with ProgressBar():
    #da.to_zarr(darr,snakemake.output.zarr,overwrite=True,dimension_separator='/', zarr_format=2)
    da.to_zarr(darr,snakemake.output.zarr,overwrite=True, zarr_format=3)


