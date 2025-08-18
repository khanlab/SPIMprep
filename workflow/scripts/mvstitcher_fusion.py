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

    print({'z':-metadata['lookup_tile_offset_z'][key],
                       'y':-metadata['lookup_tile_offset_y'][key],
                       'x':metadata['lookup_tile_offset_x'][key]})
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
#    print(msim)
    affine_from_metadata = msi_utils.get_transform_from_msim(msim, transform_key='affine_metadata')
#    print('affine from metadata')
#    print(affine_from_metadata)
    print('affine from registration')
    print(affine)

    

#output_spacing=si_utils.get_spacing_from_sim(msi_utils.get_sim_from_msim(msims[0]))

#output_spacing = {'z': 16.0, 'y': 16, 'x': 16}

#get spatial images
sims = [msi_utils.get_sim_from_msim(msim) for msim in msims]


#use affine_metadata (ie affines without registration) to get overall bounding box
params = [
    si_utils.get_affine_from_sim(sim, transform_key='affine_metadata')
    for sim in sims
]


output_spacing = si_utils.get_spacing_from_sim(sims[0])
output_stack_mode='union'

output_stack_properties = fusion.calc_fusion_stack_properties(
    sims,
    params=params,
    spacing=output_spacing,
    mode=output_stack_mode,
)


print('output_stack_properties')
print(output_stack_properties)

# Now run fusion 
fused = fusion.fuse(
    sims,
    transform_key='affine_registered',
    output_stack_properties=output_stack_properties,
    output_chunksize=256,
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



