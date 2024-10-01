import tifffile
import json
import dask.array as da
from dask.delayed import delayed
import dask.array.image 
from itertools import product
from dask.diagnostics import ProgressBar
from lib.dask_image import imread_pages
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


def get_tiff_num_pages(fs,gcs_uri):
    with fs.open(gcs_uri, 'rb') as file:
        return len(tifffile.TiffFile(file).pages)

def read_tiff_slice(fs,gcs_uri, key=0):
    """Read a single TIFF slice from GCS."""
    with fs.open(gcs_uri, 'rb') as file:
        return tifffile.imread(file, key=key)

def read_page_as_numpy(tif_file_uri, page, fs):
    """Gets a single page (i.e., 2D image) from a tif file z-stack stored in a cloud URI."""
    
    # Open the file from the cloud storage using fsspec
    with fs.open(tif_file_uri, 'rb') as f:
        # Read the file content into a buffer
        file_buffer = f.read()
    
    # Use pyvips to read from the buffer
    return pyvips.Image.new_from_buffer(file_buffer, "", page=page).numpy()



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

def build_zstack_from_single(gcs_uri,zstack_metadata,fs):
    """Build a z-stack from a single GCS URI  """

    pages = [i for i in range(zstack_metadata['n_z'])]
    lazy_arrays = [
        dask.delayed(read_tiff_slice)(fs,gcs_uri,page) for page in pages
    ]

    # Convert the list of delayed objects into a Dask array
    return da.stack([da.from_delayed(lazy_array, shape=zstack_metadata['shape'], dtype=zstack_metadata['dtype']) for lazy_array in lazy_arrays], axis=0)



#read metadata json
with open(snakemake.input.metadata_json) as fp:
    metadata = json.load(fp)


is_zstack = metadata['is_zstack']
is_tiled = metadata['is_tiled']

if is_tiled:
    if is_zstack:
        in_tif_pattern = snakemake.config["import_blaze"]["raw_tif_pattern_zstack"]
    else:
        in_tif_pattern = snakemake.config["import_blaze"]["raw_tif_pattern"]
else:
    in_tif_pattern = snakemake.config["import_blaze"]["raw_tif_pattern_notile"]




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
                    pages.append(da.from_delayed(delayed(read_page_as_numpy)('gcs://'+tif_file,i_z,fs),shape=(size_y,size_x),dtype='uint16'))
                
                zstacks.append(da.stack(pages))

            else:
                zstacks.append(build_zstack(fs.glob('gcs://'+in_tif_glob.format(tilex=tilex,tiley=tiley,prefix=metadata['prefixes'][0],channel=channel,zslice='*')),fs=fs))
            

        #have list of zstack dask arrays for the tile, one for each channel
        #stack them up and append to list of tiles
        tiles.append(da.stack(zstacks))


    #now we have list of tiles, each a dask array
    #stack them up to get our final array
    darr = da.stack(tiles)

else:
    print("single tile data not supported for tif_to_zarr_gcs")


#rescale intensities, and recast 
darr = darr * snakemake.params.intensity_rescaling
darr = darr.astype('uint16')

#now we can do the computation itself, storing to zarr
print('writing images to zarr with dask')
with ProgressBar():
    da.to_zarr(darr,snakemake.output.zarr,overwrite=True,dimension_separator='/')


