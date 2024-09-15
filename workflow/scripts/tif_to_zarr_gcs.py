import tifffile
import json
import dask.array as da
import dask.array.image 
from itertools import product
from dask.diagnostics import ProgressBar
import gcsfs

gcsfs_opts={'project': snakemake.params.storage_provider_settings['gcs'].get_settings().project,
                        'token': snakemake.input.creds}
fs = gcsfs.GCSFileSystem(**gcsfs_opts)


def replace_square_brackets(pattern):
    """replace all [ and ] in the string (have to use 
    intermediate variable to avoid conflicts)"""
    pattern = pattern.replace('[','##LEFTBRACKET##')
    pattern = pattern.replace(']','##RIGHTBRACKET##')
    pattern = pattern.replace('##LEFTBRACKET##','[[]')
    pattern = pattern.replace('##RIGHTBRACKET##','[]]')
    return pattern

def read_tiff_slice(fs,gcs_uri, key=0):
    """Read a single TIFF slice from GCS."""
    with fs.open(gcs_uri, 'rb') as file:
        return tifffile.imread(file, key=key)

def build_zstack(gcs_uris,fs):
    """Build a z-stack from a list of GCS URIs."""
    lazy_arrays = [
        dask.delayed(read_tiff_slice)(fs,uri) for uri in gcs_uris
    ]
    sample_array = read_tiff_slice(fs,gcs_uris[0])  # Read a sample to get shape and dtype
    shape = (len(gcs_uris),) + sample_array.shape
    dtype = sample_array.dtype

    # Convert the list of delayed objects into a Dask array
    return da.stack([da.from_delayed(lazy_array, shape=sample_array.shape, dtype=dtype) for lazy_array in lazy_arrays], axis=0)



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

            
        zstacks.append(build_zstack(fs.glob('gcs://'+in_tif_glob.format(tilex=tilex,tiley=tiley,prefix=metadata['prefixes'][0],channel=channel,zslice='*')),fs=fs))
        

    #have list of zstack dask arrays for the tile, one for each channel
    #stack them up and append to list of tiles
    tiles.append(da.stack(zstacks))


#now we have list of tiles, each a dask array
#stack them up to get our final array
darr = da.stack(tiles)

#rescale intensities, and recast 
darr = darr * snakemake.params.intensity_rescaling
darr = darr.astype('int16')

#now we can do the computation itself, storing to zarr
print('writing images to zarr with dask')
with ProgressBar():
    da.to_zarr(darr,snakemake.output.zarr,overwrite=True,dimension_separator='/')


