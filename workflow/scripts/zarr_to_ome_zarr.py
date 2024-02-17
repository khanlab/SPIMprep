import json
import zarr
import dask.array as da
from ome_zarr.writer import write_image
from ome_zarr.format import format_from_version
from ome_zarr.scale import Scaler
from dask.diagnostics import ProgressBar

in_zarr=snakemake.input.zarr
metadata_json=snakemake.input.metadata_json
downsampling=snakemake.params.downsampling
max_layer=snakemake.params.max_downsampling_layers #number of downsamplings by 2 to include in zarr
rechunk_size=snakemake.params.rechunk_size
out_zarr=snakemake.output.zarr
scaling_method=snakemake.params.scaling_method

#open zarr to get group name
zi = zarr.open(in_zarr)
group_name = [g for g in zi.group_keys()][0]

darr = da.from_zarr(in_zarr,component=f'{group_name}/s0',chunks=rechunk_size)
darr


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
    
    coordinate_transformations.append( [{'scale': [voxdim[0],(2**l)*voxdim[1],(2**l)*voxdim[2]], #image-pyramids in XY only
                                            'type': 'scale'}]) 


axes =  [  {'name': ax, 'type': 'space', 'unit': 'micrometer'} for ax in ['z','y','x'] ] 



store = zarr.DirectoryStore(out_zarr)
root = zarr.group(store,path='/',overwrite=True)
scaler = Scaler(max_layer=max_layer,method=scaling_method)



with ProgressBar():
    write_image(image=darr,
                            group=root,
                            scaler=scaler,
                            coordinate_transformations=coordinate_transformations,
                            storage_options={'dimension_separator': '/'},
                            axes=axes)


