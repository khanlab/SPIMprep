import json
import os
import numpy as np
from dask.diagnostics import ProgressBar
import dask.array as da
from multiview_stitcher import spatial_image_utils as si_utils
from multiview_stitcher import io, ngff_utils, msi_utils, vis_utils, registration, fusion
import matplotlib.pyplot as plt
import matplotlib
import zarr

matplotlib.use('agg')

in_zarr = snakemake.input.zarr

darr = da.from_zarr(in_zarr)

(n_tiles,n_chans,n_z,n_x,n_y) = darr.shape

#read metadata json
with open(snakemake.input.metadata_json) as fp:
    metadata = json.load(fp)


msims = []
zarr_paths = []


channel=metadata['channels'][snakemake.params.channel_index]

for i_tile in metadata['chunks']:
 
    key = f'tile-{i_tile}_chan-{channel}_z-0000'

#    zarr_path =f'tile_{i_tile:02d}.zarr'

    # read tile image
    im_data = da.squeeze(darr[i_tile,:,:,:,:]) 
       
    print(im_data.shape)
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

    # write to OME-Zarr
 #   print(f'writing tile {i_tile} to ome_zarr') 
 #   with ProgressBar():
 #       ngff_utils.write_sim_to_ome_zarr(sim, zarr_path)


    #for next steps, we read things back as msim
#    sim = ngff_utils.read_sim_from_ome_zarr(zarr_path)
    msim = msi_utils.get_msim_from_sim(sim)
#    zarr_paths.append(zarr_path)

    msims.append(msim)


"""
print('saving tile visualization')
fig, ax = vis_utils.plot_positions(
    msims,
    use_positional_colors=True, # set to False for faster execution in case of more than 20 tiles/views
    transform_key='affine_metadata'
    )


plt.savefig(snakemake.output.tiling_qc_png)
"""

# interrupt snakemake to stop the viewer
#vis_utils.view_neuroglancer(
#    sims=[msi_utils.get_sim_from_msim(msim) for msim in msims],
#    ome_zarr_paths=zarr_paths,
#    transform_key=io.METADATA_TRANSFORM_KEY,
#    port=8001
#)

curr_transform_key = 'affine_metadata'

print('performing stitching registration')
with ProgressBar():
    params = registration.register(
        msims,
        reg_channel_index=snakemake.params.reg_channel_index,
        transform_key=curr_transform_key,
        new_transform_key='affine_registered',
        pre_registration_pruning_method="keep_axis_aligned", # works well for tiles on a grid
        scheduler="threads",
        **snakemake.params.registration_opts,
    )


for imsim, msim in enumerate(msims):
    affine = np.array(msi_utils.get_transform_from_msim(msim, transform_key='affine_registered')[0])
    print(f'tile index {imsim}\n', affine)


fused = fusion.fuse(
    [msi_utils.get_sim_from_msim(msim) for msim in msims],
    transform_key='affine_registered',
    fusion_func=fusion.max_fusion,
    #    **snakemake.params.fusion_opts,
    )

print(fused)
print(type(fused))
print(fused.shape)

print('shape of array to save')
print(fused.data[0].shape)

#save each channel separately (to be consistent with legacy bigstitcher workflow)
with ProgressBar():
    fused.data[0][snakemake.params.channel_index,:,:,:].to_zarr(snakemake.output.zarr,overwrite=True,
                                                                    dimension_separator='/',component='fused/s0',zarr_format=2)

"""
znimg = ZarrNii.from_darr(fused.data[0], axes_order='ZYX',spacing=( metadata['physical_size_z'],
               metadata['physical_size_y'],
               metadata['physical_size_x']))




print(f'Fusing views and saving output to ome zarr...')
with ProgressBar():
    fused = ngff_utils.write_sim_to_ome_zarr(
        fused, snakemake.output.zarr, overwrite=True
    )
"""
