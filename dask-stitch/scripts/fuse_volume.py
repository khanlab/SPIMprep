from lib.utils import *


# Paths to input tiles
ome_zarr_paths = snakemake.input.ome_zarr

# Build R-tree index for all chunks in all tiles
rtree_idx, znimgs = build_rtree_index(ome_zarr_paths)


# Parameters for the fused array
fused_shape = (128,128,128)  # Full shape of the fused array
chunk_shape = (64,64,64)  # Chunk size of the fused array
global_bbox_min = np.array([0, 0, 0])  # Minimum physical coordinate of the fused array
target_voxel_spacing = np.array([1.0, 1.0, 1.0])  # Voxel spacing of the fused array

# Create an empty Dask array for the fused output
fused_dask_array = da.zeros(fused_shape, chunks=chunk_shape, dtype=np.float32)



# Apply map_blocks to process the fused array
processed_dask_array = fused_dask_array.map_blocks(
    process_block,
    dtype=np.float32,
    rtree_idx=rtree_idx,
    znimgs=znimgs,
    global_bbox_min=global_bbox_min,
    target_voxel_spacing=target_voxel_spacing,
)


znimg_out = ZarrNii.from_darr(processed_dask_array.reshape((1,processed_dask_array.shape[0],processed_dask_array.shape[1],processed_dask_array.shape[2])))
znimg_out.to_ome_zarr(snakemake.output.ome_zarr)
znimg_out.to_nifti(snakemake.output.nifti)

"""

print('doing test query')
# Define the bounding box of an output chunk (example)
output_chunk_bbox =( -1, -1, -1, 0,0,0)  # Replace with actual bounding box of a fused chunk
output_chunk_bbox =( -36, -33, -49, -35,-32,-47)  # Replace with actual bounding box of a fused chunk
print(output_chunk_bbox)

# Query intersecting chunks
intersecting_chunks = query_intersections(rtree_idx, output_chunk_bbox)

print(f'found {len(intersecting_chunks)} chunks')
# Extract the attached data (tile_idx, chunk_id)
for entry in intersecting_chunks:
    tile_idx, chunk_id = entry.object
    print(f"Intersecting chunk: tile_idx={tile_idx}, chunk_id={chunk_id}")

"""
