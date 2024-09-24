import tifffile
import json
import dask.array as da
import dask.array.image 
from itertools import product
from dask.diagnostics import ProgressBar

def replace_square_brackets(pattern):
    """replace all [ and ] in the string (have to use 
    intermediate variable to avoid conflicts)"""
    pattern = pattern.replace('[','##LEFTBRACKET##')
    pattern = pattern.replace(']','##RIGHTBRACKET##')
    pattern = pattern.replace('##LEFTBRACKET##','[[]')
    pattern = pattern.replace('##RIGHTBRACKET##','[]]')
    return pattern


def single_imread(*args):
    """create function handle to tifffile.imread 
    that sets key=0"""
    return tifffile.imread(*args,key=0)


#use tif pattern but replace the [ and ] with [[] and []] so glob doesn't choke
in_tif_glob = replace_square_brackets(str(snakemake.params.in_tif_pattern))


#read metadata json
with open(snakemake.input.metadata_json) as fp:
    metadata = json.load(fp)

#TODO: put these in top-level metadata for easier access..
size_x=metadata['ome_full_metadata']['OME']['Image']['Pixels']['@SizeX']
size_y=metadata['ome_full_metadata']['OME']['Image']['Pixels']['@SizeY']
size_z=metadata['ome_full_metadata']['OME']['Image']['Pixels']['@SizeZ']
size_c=metadata['ome_full_metadata']['OME']['Image']['Pixels']['@SizeC']
size_tiles=len(metadata['tiles_x'])*len(metadata['tiles_y'])


#now get the first channel and first zslice tif 
tiles=[]
for i_tile,(tilex,tiley) in enumerate(product(metadata['tiles_x'],metadata['tiles_y'])):
        
    zstacks=[]
    for i_chan,channel in enumerate(metadata['channels']):

            
        zstacks.append(dask.array.image.imread(in_tif_glob.format(tilex=tilex,tiley=tiley,prefix=metadata['prefixes'][0],channel=channel,zslice='*'), imread=single_imread))
        

    #have list of zstack dask arrays for the tile, one for each channel
    #stack them up and append to list of tiles
    tiles.append(da.stack(zstacks))


#now we have list of tiles, each a dask array
#stack them up to get our final array
darr = da.stack(tiles)

#rescale intensities, and recast 
darr = darr * snakemake.params.intensity_rescaling
darr = darr.astype('uint16')

#now we can do the computation itself, storing to zarr
print('writing images to zarr with dask')
with ProgressBar():
    da.to_zarr(darr,snakemake.output.zarr,overwrite=True,dimension_separator='/')


