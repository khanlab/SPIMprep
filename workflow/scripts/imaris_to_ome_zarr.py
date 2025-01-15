import h5py
import hdf5plugin
import zarr
import json
import zarr
import dask.array as da
from ome_zarr.writer import write_image
from ome_zarr.format import format_from_version
from ome_zarr.scale import Scaler
from dask.diagnostics import ProgressBar
from upath import UPath as Path
from lib.cloud_io import get_fsspec, is_remote



stains=snakemake.params.stains
metadata_json=snakemake.input.metadata_json
downsampling=snakemake.params.downsampling
max_layer=snakemake.params.max_downsampling_layers #number of downsamplings by 2 to include in zarr
rechunk_size=snakemake.params.rechunk_size
out_zarr=snakemake.output.zarr
scaling_method=snakemake.params.scaling_method
uri = snakemake.params.uri



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


if is_remote(uri):
    fs_args={'storage_provider_settings':snakemake.params.storage_provider_settings,'creds':snakemake.input.creds}
    fs = get_fsspec(uri,**fs_args)
    store = zarr.storage.FSStore(Path(uri).path,fs=fs,dimension_separator='/',mode='w')
else:
    #store = zarr.DirectoryStore(out_zarr,dimension_separator='/') 
    store = zarr.ZipStore(out_zarr,dimension_separator='/',mode='x') 



with h5py.File(snakemake.input.ims, 'r') as h5_file:

    darr_list=[]
    for zarr_i,stain in enumerate(stains):
        h5_group=f'DataSet/ResolutionLevel 0/TimePoint 0/Channel {zarr_i}/Data'
        darr_list.append(da.from_array(h5_file[h5_group]).rechunk(rechunk_size))


        #append to omero metadata
        channel_metadata={key:val for key,val in snakemake.config['ome_zarr']['omero_metadata']['channels']['defaults'].items()}
        channel_name=stain
        channel_metadata['label'] = channel_name
        default_color=snakemake.config['ome_zarr']['omero_metadata']['channels']['default_color']
        color=snakemake.config['ome_zarr']['omero_metadata']['channels']['color_mapping'].get(channel_name,default_color)
        channel_metadata['color'] = color
        omero['channels'].append(channel_metadata)
        

    darr_channels = da.stack(darr_list)



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


