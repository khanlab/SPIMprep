import json
import zarr
import dask.array as da
from dask.array.image import imread as dask_imread
from ome_zarr.io import parse_url
from ome_zarr.writer import write_image
from ome_zarr.format import format_from_version
from ome_zarr.scale import Scaler
from dask.diagnostics import ProgressBar

in_tif_glob = snakemake.params.in_tif_glob

metadata_json=snakemake.input.metadata_json
downsampling=snakemake.params.downsampling
max_layer=snakemake.params.max_downsampling_layers #number of downsamplings by 2 to include in zarr
rechunk_size=snakemake.params.rechunk_size
out_zarr=snakemake.output.zarr
stains=snakemake.params.stains
scaling_method=snakemake.params.scaling_method

# prepare metadata for ome-zarr
with open(metadata_json) as fp:
    metadata = json.load(fp)

voxdim = [float(metadata['physical_size_z']) * float(downsampling) / 1000.0,
            float(metadata['physical_size_y']) * float(downsampling) / 1000.0,
            float(metadata['physical_size_x']) * float(downsampling) / 1000.0]


coordinate_transformations = []
#for each resolution (dataset), we have a list of dicts, transformations to apply.. 
#in this case just a single one (scaling by voxel size)

for l in range(max_layer+1):
    
    coordinate_transformations.append( [{'scale': [1,voxdim[0],(2**l)*voxdim[1],(2**l)*voxdim[2]], #image-pyramids in XY only
                                            'type': 'scale'}]) 


axes =  [{'name': 'c', 'type': 'channel'}] + [{'name': ax, 'type': 'space', 'unit': 'millimeter'} for ax in ['z','y','x'] ] 


#init omero metadata
omero={key:val for key,val in snakemake.config['ome_zarr']['omero_metadata']['defaults'].items()}
omero['channels']=[]

darr_list=[]
for i,stain in enumerate(stains):
    
    tif_glob = in_tif_glob.format(stain=stain)
    darr_list.append(dask_imread(tif_glob).rechunk(rechunk_size))

    #append to omero metadata
    channel_metadata={key:val for key,val in snakemake.config['ome_zarr']['omero_metadata']['channels']['defaults'].items()}
    channel_name=stains[i]
    channel_metadata['label'] = channel_name
    default_color=snakemake.config['ome_zarr']['omero_metadata']['channels']['default_color']
    color=snakemake.config['ome_zarr']['omero_metadata']['channels']['color_mapping'].get(channel_name,default_color)
    channel_metadata['color'] = color
    omero['channels'].append(channel_metadata)
    

darr_channels = da.stack(darr_list)


if out_zarr.endswith('zarr.touch'):
    #use the uri
    uri = snakemake.params.uri
    print(uri)

    #strip off gcs:// for use in gcsfs
    if uri.startswith('gcs://'):
        uri = uri[6:]
        import gcsfs
        gcsfs_opts={'project': 't-system-193821'}
        fs = gcsfs.GCSFileSystem(**gcsfs_opts)
        store = zarr.storage.FSStore(uri,fs=fs,dimension_separator='/',mode='w')
    else:
        print(f'cannot parse uri {uri}')
    zi = zarr.open(store=store,mode='w')


   
else:
    #load as directorystore
    store = zarr.DirectoryStore(out_zarr)
    zi = zarr.open(out_zarr,mode='w')



group = zarr.group(store,overwrite=True)
scaler = Scaler(max_layer=max_layer,method=scaling_method)


with ProgressBar():
    write_image(image=darr_channels,
                            group=group,
                            scaler=scaler,
                            coordinate_transformations=coordinate_transformations,
                            axes=axes,
                            metadata={'omero':omero}
                                )


