from lib.utils import *
from zarrnii import ZarrNii

# Paths to input tiles
ome_zarr_paths = snakemake.input.ome_zarr

# Build R-tree index for all chunks in all tiles
rtree_idx, znimgs = build_rtree_index(ome_zarr_paths)


# Parameters for the fused array
output_shape=(400,400,400)
output_chunks=(100,100,100)
fused_dask_array = da.zeros(output_shape, chunks=output_chunks, dtype=np.float32)
znimg_out = ZarrNii.from_darr(fused_dask_array,
                                             axes_order=znimgs[0].axes_order,
                                             orientation=znimgs[0].get_orientation(),
                                             spacing=(1,1,1),
                                             origin=(-100,-100,-100))

znimg_fused = fuse_output_blocks(znimg_out,rtree_idx=rtree_idx,tile_znimgs=znimgs)

znimg_fused.to_ome_zarr(snakemake.output.ome_zarr)

ZarrNii.from_ome_zarr(snakemake.output.ome_zarr).to_nifti(snakemake.output.nifti)



