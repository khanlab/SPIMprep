import zarr
from dask.diagnostics import ProgressBar
from upath import UPath as Path
from zarrnii import ZarrNii

uri = snakemake.params.uri
in_zarr = snakemake.input.zarr
channel_index = snakemake.params.channel_index


if Path(uri).suffix == '.zip':
    store = zarr.storage.ZipStore(Path(uri).path,mode='r')
else:
    store = zarr.storage.LocalStore(Path(uri).path,read_only=True)


znimg = ZarrNii.from_ome_zarr(
    store, level=int(snakemake.wildcards.level), channels=[channel_index]
)


# before updating zarrnii ngffzarr3 branch to accommodate anisotropically downsampled data, instead
# we will calculate the z downsampling factor and downsample accordingly - TODO: move this to zarrnii

import numpy as np

# Get scale and axes order
scale = znimg.coordinate_transformations[0].scale

axes = znimg.axes  # list of Axis objects

# Build a mapping from axis name to index
axis_index = {axis.name.lower(): i for i, axis in enumerate(axes)}

# Extract x and z scales
x_scale = scale[axis_index["x"]]
z_scale = scale[axis_index["z"]]

# Compute ratio and power
ratio = x_scale / z_scale
level = int(np.log2(round(ratio)))


with ProgressBar():
    if level == 0:
        znimg.to_nifti(snakemake.output.nii)
    else:
        znimg.downsample(along_z=2**level).to_nifti(snakemake.output.nii)


"""

if is_remote(uri):
    fs_args={'storage_provider_settings':snakemake.params.storage_provider_settings,'creds':snakemake.input.creds}
else:
    fs_args={}

fs = get_fsspec(uri,**fs_args)

if Path(uri).suffix == '.zip':
    store = zarr.storage.ZipStore(Path(uri).path,dimension_separator='/',mode='r')
else:
    store = zarr.storage.FSStore(Path(uri).path,fs=fs,dimension_separator='/',mode='r')

zi = zarr.open(store=store,mode='r')
 

attrs=zi['/'].attrs.asdict()

level=int(snakemake.wildcards.level)

#read coordinate transform from ome-zarr
transforms = attrs['multiscales'][0]['datasets'][level]['coordinateTransformations']


#zarr uses z,y,x ordering, we reverse this for nifti
# also flip to set orientation properly
affine = np.eye(4)
affine[0,0]=-transforms[0]['scale'][3] #x
affine[1,1]=-transforms[0]['scale'][2] #y
affine[2,2]=-transforms[0]['scale'][1] #z

if is_remote(uri):
    #grab the channel index corresponding to the stain
    darr = da.from_zarr(store,component=f'/{level}')[channel_index,:,:,:].squeeze()
else:
    darr = da.from_zarr(in_zarr,component=f'/{level}')[channel_index,:,:,:].squeeze()

#input array axes are ZYX 
#writing to nifti we want XYZ
out_arr = np.moveaxis(darr,(0,1,2),(2,1,0)) 

nii = nib.Nifti1Image(out_arr,
                    affine=affine
                    )
                    
nii.to_filename(snakemake.output.nii)
"""
