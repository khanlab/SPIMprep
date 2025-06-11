import json
import os
import numpy as np
from dask.diagnostics import ProgressBar
import dask.array as da
from multiview_stitcher import spatial_image_utils as si_utils
from multiview_stitcher import io, ngff_utils, msi_utils, vis_utils, registration, fusion, param_utils
import matplotlib.pyplot as plt
import matplotlib
import zarr

matplotlib.use('agg')
import dask

dask.config.set(scheduler='threads', num_workers=min(32,snakemake.threads))


in_zarr = snakemake.input.zarr

darr = da.from_zarr(in_zarr)

(n_tiles,n_chans,n_z,n_x,n_y) = darr.shape

#read metadata json
with open(snakemake.input.metadata_json) as fp:
    metadata = json.load(fp)


msims = []
zarr_paths = []


reg_channel=metadata['channels'][snakemake.params.reg_channel_index]

for i_tile in metadata['chunks']:
 
    key = f'tile-{i_tile}_chan-{reg_channel}_z-0000'

    # read tile image
    im_data = da.squeeze(darr[i_tile,:,:,:,:]) 


    print(f'shape: {im_data.shape}, chunksize: {im_data.chunksize}')
    sim = si_utils.get_sim_from_array(
        im_data,
        dims=["c","z","y","x"],
        scale={'z': metadata['physical_size_z'],
               'y':metadata['physical_size_y'],
               'x':metadata['physical_size_x']},
        translation = {'z':-metadata['lookup_tile_offset_z'][key],
                       'y':-metadata['lookup_tile_offset_y'][key],
                       'x':metadata['lookup_tile_offset_x'][key]},
        c_coords=snakemake.params.channels,
        transform_key=io.METADATA_TRANSFORM_KEY,
        )

    #for next steps, we read things back as msim
    msim = msi_utils.get_msim_from_sim(sim)
    msims.append(msim)



# Load affines from disk
affines_npz = np.load(snakemake.input.affines)
# Example: affines_npz['tile_0'], affines_npz['tile_1'], ...

# Restore affines into msims
for imsim, msim in enumerate(msims):
    affine = affines_npz[f"tile_{imsim}"]
    # Set the affine back in the msim metadata
    affine_xr = param_utils.affine_to_xaffine(affine)
    msi_utils.set_affine_transform(msim, affine_xr, transform_key='affine_registered')

#output_spacing=si_utils.get_spacing_from_sim(msi_utils.get_sim_from_msim(msims[0]))

#output_spacing = {'z': 16.0, 'y': 16, 'x': 16}

# Now run fusion as before
fused = fusion.fuse(
    [msi_utils.get_sim_from_msim(msim) for msim in msims],
    transform_key='affine_registered',
#    output_stack_mode='sample',
#    fusion_func=fusion.max_fusion,
#    output_spacing=output_spacing,
    output_chunksize=256,
    #    **snakemake.params.fusion_opts,
)

print(fused)
print(type(fused))
print(fused.shape)

print('shape of array to save')
print(fused.data[0].shape)




with ProgressBar():
    fused = ngff_utils.write_sim_to_ome_zarr(
        fused, snakemake.output.zarr, overwrite=True
    )



##save individual channel (to be consistent with legacy bigstitcher workflow)
#with ProgressBar():
#    fused.data[0][snakemake.params.channel_index,:,:,:].to_zarr(snakemake.output.zarr,overwrite=True,
#                                                                    dimension_separator='/',component='fused/s0',zarr_format=2)

