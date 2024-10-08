import tifffile
import os
import json
import pyvips
import dask.array as da
from dask.delayed import delayed
from dask.array.image import imread as imread_tifs
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


def read_page_as_numpy(tif_file,page):    
    """gets a single page (i.e. 2d image) from a tif file zstack"""
    return pyvips.Image.new_from_file(tif_file, page=page).numpy()


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


#TODO: put these in top-level metadata for easier access..
size_x=int(metadata['ome_full_metadata']['OME']['Image']['Pixels']['@SizeX'])
size_y=int(metadata['ome_full_metadata']['OME']['Image']['Pixels']['@SizeY'])
size_z=int(metadata['ome_full_metadata']['OME']['Image']['Pixels']['@SizeZ'])
size_c=int(metadata['ome_full_metadata']['OME']['Image']['Pixels']['@SizeC'])




if is_tiled:

    size_tiles=len(metadata['tiles_x'])*len(metadata['tiles_y'])

    #now get the first channel and first zslice tif 
    tiles=[]
    for i_tile,(tilex,tiley) in enumerate(product(metadata['tiles_x'],metadata['tiles_y'])):
        print(f'tile {tilex}x{tiley}, {i_tile}') 
        zstacks=[]
        for i_chan,channel in enumerate(metadata['channels']):

            print(f'channel {i_chan}')

            if is_zstack:
                tif_file = in_tif_pattern.format(tilex=tilex,tiley=tiley,prefix=metadata['prefixes'][0],channel=channel)
                
                pages=[]
                #read each page
                for i_z in range(size_z):
                    pages.append(da.from_delayed(delayed(read_page_as_numpy)(tif_file,i_z),shape=(size_y,size_x),dtype='uint16'))
                
                zstacks.append(da.stack(pages))
                print(zstacks[-1].shape)
            else:
                zstacks.append(imread_tifs(in_tif_glob.format(tilex=tilex,tiley=tiley,prefix=metadata['prefixes'][0],channel=channel,zslice='*'), imread=single_imread))

        #have list of zstack dask arrays for the tile, one for each channel
        #stack them up and append to list of tiles
        tiles.append(da.stack(zstacks))



    #now we have list of tiles, each a dask array
    #stack them up to get our final array
    darr = da.stack(tiles)

else:
    #single tile, zslices:
    zstacks=[]
    for i_chan,channel in enumerate(metadata['channels']):

        print(f'channel {i_chan}')
        zstacks.append(imread_tifs(in_tif_glob.format(prefix=metadata['prefixes'][0],channel=channel,zslice='*'), imread=single_imread))

    #stack the channels up
    darr = da.stack(zstacks)



print(darr.shape)
print(darr.chunksize)

#rescale intensities, and recast 
darr = darr * snakemake.params.intensity_rescaling
darr = darr.astype('uint16')
#now we can do the computation itself, storing to zarr
print('writing images to zarr with dask')
with ProgressBar():
    da.to_zarr(darr,snakemake.output.zarr,overwrite=True,dimension_separator='/')


