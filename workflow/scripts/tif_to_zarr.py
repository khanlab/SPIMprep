import tempfile
import numpy as np
import tifffile
import os
import json
import pyvips
import zarr
import dask
import dask.array as da
from dask.delayed import delayed
from dask.array.image import imread as imread_tifs
from itertools import product
from dask.diagnostics import ProgressBar
from pathlib import Path

#this limits workers to number of threads assigned by snakemake
dask.config.set(scheduler='threads', num_workers=snakemake.threads)  


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

def read_stack_as_numpy(tif_file, Nz,Ny,Nx):
    """Gets the full stack (i.e., 3D image) from a tif file z-stack stored in a cloud URI."""
    
    #init array
    vol = np.zeros((Nz,Ny,Nx),dtype='uint16')

    # Use pyvips to read from the file
    for i in range(Nz):
        vol[i,:,:] = pyvips.Image.new_from_file(tif_file,  page=i).numpy()

    return vol


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
                
                zstacks.append(da.from_delayed(delayed(read_stack_as_numpy)(tif_file,size_z,size_y,size_x),shape=(size_z,size_y,size_x),dtype='uint16').rechunk((1,size_y,size_x)))

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
if Path(snakemake.output.zarr).suffix == '.zip':
    store = zarr.storage.ZipStore(snakemake.output.zarr,dimension_separator='/',mode='x')
else:
    store = zarr.storage.DirectoryStore(snakemake.output.zarr,dimension_separator='/')

print('writing images to zarr with dask')
with ProgressBar():
    da.to_zarr(darr,store,overwrite=True)


