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


def convert_hdf5_to_zarr(hdf5_path, zarr_path, chunks):
    """
    Convert an HDF5 file to Zarr using h5py and zarr, handling chunked copying.

    Parameters:
        hdf5_path (str): Path to the input HDF5 (.ims) file.
        zarr_path (str): Path to the output Zarr dataset.
        chunks (tuple): Chunk size for the Zarr dataset.
    """

    h5py._errors.unsilence_errors()
    # Open the HDF5 file and create a Zarr root group
    with h5py.File(hdf5_path, "r") as hdf5_file:
        zarr_store = zarr.open_group(zarr_path, mode="w")

        # Define the specific path to copy
        target_path = "DataSet/ResolutionLevel 0/TimePoint 0"

        # Check if the target path exists in HDF5
        if target_path in hdf5_file:
            hdf5_group = hdf5_file[target_path]

            def copy_group(hdf5_group, zarr_group):
                """
                Copies channel groups and their 'Data' datasets chunk by chunk.

                Args:
                    hdf5_group: HDF5 group containing the dataset.
                    zarr_group: Zarr group to write the dataset to.
                """
                for key, item in hdf5_group.items():
                    if isinstance(item, h5py.Group) and key.startswith("Channel"):  # Only copy Channel groups
                        channel_group = item
                        if "Data" in channel_group:  # Only copy the Data dataset in each Channel
                            data_item = channel_group["Data"]

                            # Create the Zarr dataset
                            zarr_dataset = zarr_group.require_dataset(
                                name=key + "/Data",
                                shape=data_item.shape,
                                chunks=chunks,
                                dtype=data_item.dtype,
                                compression="blosc",  # Optional compression
                            )

                            # Copy data chunk by chunk
                            for i_start in range(0, data_item.shape[0], chunks[0]):
                                for j_start in range(0, data_item.shape[1], chunks[1]):
                                    for k_start in range(0, data_item.shape[2], chunks[2]):
                                        i_end = min(i_start + chunks[0], data_item.shape[0])
                                        j_end = min(j_start + chunks[1], data_item.shape[1])
                                        k_end = min(k_start + chunks[2], data_item.shape[2])

                                        slices = (
                                            slice(i_start, i_end),
                                            slice(j_start, j_end),
                                            slice(k_start, k_end),
                                        )
                                        print(f"Copying slice {slices} for {key}")
                                        zarr_dataset[slices] = data_item[slices]

            # Start copying only the Channel groups
            copy_group(hdf5_group, zarr_store)



rechunk_size=snakemake.params.rechunk_size

#copy imaris (hdf5) to zarr
convert_hdf5_to_zarr(
    hdf5_path=snakemake.input.ims,
    zarr_path='copy_hdf5.zarr',
    chunks=rechunk_size
)



in_zarr='copy_hdf5.zarr'
metadata_json=snakemake.input.metadata_json
downsampling=snakemake.params.downsampling
max_layer=snakemake.params.max_downsampling_layers #number of downsamplings by 2 to include in zarr
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
    darr_list.append(da.from_zarr(in_zarr,component=f'Channel {zarr_i}/Data'))


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


