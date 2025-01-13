import h5py
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


def convert_hdf5_to_zarr(hdf5_path, zarr_path):
    """
    Convert an HDF5 file to Zarr using h5py and zarr.

    Parameters:
        hdf5_path (str): Path to the input HDF5 (.ims) file.
        zarr_path (str): Path to the output Zarr dataset.
    """
    # Open the HDF5 file and create a Zarr root group
    with h5py.File(hdf5_path, "r") as hdf5_file:
        zarr_store = zarr.open_group(zarr_path, mode="w")

        def copy_group(hdf5_group, zarr_group):
            for key, item in hdf5_group.items():
                if isinstance(item, h5py.Group):  # Recursively copy groups
                    new_group = zarr_group.create_group(key)
                    copy_group(item, new_group)
                elif isinstance(item, h5py.Dataset):  # Copy datasets
                    zarr_group.create_dataset(
                        name=key,
                        data=item[()],
                        chunks=item.chunks,
                        dtype=item.dtype,
                        compression="blosc"  # Optional compression
                    )
                    print(f"Copied dataset: {key}")

        # Start copying from the root group
        copy_group(hdf5_file, zarr_store)

    print(f"Converted HDF5 file to Zarr at: {zarr_path}")


#copy imaris (hdf5) to zarr -- TODO: don't need to copy everything 
convert_hdf5_to_zarr(
    hdf5_path=snakemake.input.ims,
    zarr_path='copy_hdf5.zarr',
)



in_zarr='copy_hdf5.zarr'
metadata_json=snakemake.input.metadata_json
downsampling=snakemake.params.downsampling
max_layer=snakemake.params.max_downsampling_layers #number of downsamplings by 2 to include in zarr
rechunk_size=snakemake.params.rechunk_size
out_zarr=snakemake.output.zarr
stains=snakemake.params.stains
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
    store = zarr.DirectoryStore(out_zarr,dimension_separator='/') 





darr_list=[]
for zarr_i,stain in enumerate(stains):
    #open zarr to get group name
    zi = zarr.open(in_zarr)
#    darr_list.append(da.from_zarr(in_zarr,component=f'DataSet/ResolutionLevel 0/TimePoint 0/Channel {zarr_i}/Data',chunks=rechunk_size))
    darr_list.append(da.from_zarr(in_zarr,component=f'DataSet/ResolutionLevel 0/TimePoint 0/Channel {zarr_i}/Data'))


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


