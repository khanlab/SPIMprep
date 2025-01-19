import numpy as np
import dask.array as da
from zarrnii import ZarrNii
from rtree import index
import dask
dask.config.set(scheduler='synchronous')  # overwrite default with single-threaded scheduler


def compute_chunk_bounding_boxes(dask_array, affine, tile_index):
    """
    Compute the bounding boxes of each chunk in a Dask array in physical space.

    Parameters:
    - dask_array (dask.array): The input tile as a Dask array.
    - affine (np.ndarray): 4x4 affine matrix mapping voxel indices to physical space.
    - tile_index (int): Index of the tile for tracking.

    Returns:
    - List of bounding boxes in physical space.
      Each bounding box is (chunk_id, (min_x, min_y, min_z, max_x, max_y, max_z)).
    """

    chunk_bboxes = []
    chunk_shapes = dask_array.chunks  # Tuple of chunk sizes along each dimension

    # Iterate through all chunk indices
    for z_idx, z_start in enumerate(np.cumsum((0,) + chunk_shapes[0][:-1])):
        z_end = z_start + chunk_shapes[0][z_idx]

        for y_idx, y_start in enumerate(np.cumsum((0,) + chunk_shapes[1][:-1])):
            y_end = y_start + chunk_shapes[1][y_idx]

            for x_idx, x_start in enumerate(np.cumsum((0,) + chunk_shapes[2][:-1])):
                x_end = x_start + chunk_shapes[2][x_idx]

                # Voxel coordinates for the chunk
                voxel_min = np.array([z_start, y_start, x_start])
                voxel_max = np.array([z_end, y_end, x_end])

                # Convert voxel coordinates to physical space using affine
                physical_corners = np.array([
                    affine @ np.array([voxel_min[2], voxel_min[1], voxel_min[0], 1]),
                    affine @ np.array([voxel_max[2], voxel_min[1], voxel_min[0], 1]),
                    affine @ np.array([voxel_min[2], voxel_max[1], voxel_min[0], 1]),
                    affine @ np.array([voxel_min[2], voxel_min[1], voxel_max[0], 1]),
                    affine @ np.array([voxel_max[2], voxel_max[1], voxel_min[0], 1]),
                    affine @ np.array([voxel_min[2], voxel_max[1], voxel_max[0], 1]),
                    affine @ np.array([voxel_max[2], voxel_min[1], voxel_max[0], 1]),
                    affine @ np.array([voxel_max[2], voxel_max[1], voxel_max[0], 1]),
                ])

                physical_bbox_min = physical_corners.min(axis=0)[:3]
                physical_bbox_max = physical_corners.max(axis=0)[:3]

                # Store the bounding box with a unique chunk_id
                chunk_id = (z_idx, y_idx, x_idx)
                chunk_bboxes.append((chunk_id, (physical_bbox_min, physical_bbox_max, tile_index)))

    return chunk_bboxes


def build_rtree_index(ome_zarr_paths):
    """
    Build an R-tree index of all chunk bounding boxes for all tiles.

    Parameters:
    - ome_zarr_paths (list of str): Paths to all input tiles as OME-Zarr datasets.

    Returns:
    - rtree.index.Index: R-tree index for fast spatial queries.
    - List of ZarrNii objects corresponding to the input paths.
    """
    idx = index.Index()
    p = index.Property()
    p.dimension = 3
#    p.dat_extension = 'data' #this is to persist to disk
#    p.idx_extension = 'index'
    #idx3d = index.Index('3d_index',properties=p)
    idx3d = index.Index(properties=p)

    znimgs = []  # Store the loaded ZarrNii objects for later use

    for tile_index, path in enumerate(ome_zarr_paths):
        #print(path)
        #print(tile_index)
        znimg = ZarrNii.from_ome_zarr(path)
        znimgs.append(znimg)

        dask_array = znimg.darr.squeeze()  # Assume the data is squeezed to 3D
        affine = znimg.affine

        #HACK fix
#        affine[:3,3] = -1 * np.flip(affine[:3,3])


#        print('blocks')
#        print(dask_array.blocks)
#        print('chunks')
#        print(dask_array.chunks)

        # Compute bounding boxes for each chunk
        chunk_bboxes = compute_chunk_bounding_boxes(dask_array, affine, tile_index)
#        print('first bbox')
#        print(chunk_bboxes[0])
#        print('last bbox')
#        print(chunk_bboxes[-1])
#        print('all ')
#        print(chunk_bboxes)
        # Insert bounding boxes into the R-tree
        for chunk_id, (physical_bbox_min, physical_bbox_max, tile_idx) in chunk_bboxes:
            print(f'{chunk_id}, {physical_bbox_min} to {physical_bbox_max} in tile {tile_idx}')
            rtree_bbox = (
                physical_bbox_min[0], physical_bbox_min[1], physical_bbox_min[2],
                physical_bbox_max[0], physical_bbox_max[1], physical_bbox_max[2],
            )
            unique_chunk_id = hash(f'{chunk_id}_{tile_idx}')  # Use hash to create a unique ID
            idx3d.insert(unique_chunk_id, rtree_bbox, obj=(tile_idx, chunk_id))


    return idx3d, znimgs


def query_intersections(rtree_idx, chunk_bbox):
    """
    Query the R-tree index for chunks that intersect with a given bounding box.

    Parameters:
    - rtree_idx (rtree.index.Index): The R-tree index.
    - chunk_bbox (tuple): The bounding box of the chunk to query (min_x, min_y, min_z, max_x, max_y, max_z).

    Returns:
    - List of intersecting chunk entries from the R-tree.
    """
    return list(rtree_idx.intersection(chunk_bbox, objects=True))



import numpy as np
from scipy.ndimage import map_coordinates



# Define the process_block function
def process_block(block_data, block_info=None, rtree_idx=None, znimgs=None, global_bbox_min=None, target_voxel_spacing=None):
    """
    Process a block of the fused array by loading and resampling overlapping chunks.

    Parameters:
    - block_data (np.ndarray): Placeholder array for the output block (filled with zeros).
    - block_info (dict): Information about the block being processed (from Dask).
    - rtree_idx (rtree.index.Index): R-tree index of the tile chunks.
    - znimgs (list): List of ZarrNii objects for the tiles.
    - global_bbox_min (np.ndarray): Minimum physical coordinates of the fused array.
    - target_voxel_spacing (np.ndarray): Target voxel spacing for the fused array.

    Returns:
    - np.ndarray: The processed block of the fused array.
    """
    # Extract block location and shape
    block_shape = block_data.shape
    block_location = block_info[0]["chunk-location"]  # (z, y, x) chunk indices
    block_bbox_min = global_bbox_min + block_location * target_voxel_spacing
    block_bbox_max = block_bbox_min + np.array(block_shape) * target_voxel_spacing
    block_bbox = (*block_bbox_min, *block_bbox_max)

    print(f'at block {block_location}')
    print(f'bbox: {block_bbox}')


    # Query the R-tree for overlapping chunks
    intersecting_chunks = list(rtree_idx.intersection(block_bbox, objects=True))

    # Initialize output block and weight array
    fused_block = np.zeros(block_shape, dtype=np.float32)
    weight = np.zeros(block_shape, dtype=np.float32)

    # Process each overlapping chunk
    for entry in intersecting_chunks:
        #print("have an intersecting chunk")
        tile_idx, chunk_id = entry.object
        znimg = znimgs[tile_idx]
        tile_dask_array = znimg.darr.squeeze()
        affine = znimg.vox2ras.affine
    
        #HACK fix
#       affine[:3,3] = -1 * np.flip(affine[:3,3])


        # Extract chunk data
        z_idx, y_idx, x_idx = chunk_id
        tile_chunk = tile_dask_array.blocks[z_idx, y_idx, x_idx].compute()

        # Compute the bounding box of the chunk in physical space
        chunk_shape = tile_chunk.shape
        voxel_min = np.array([z_idx, y_idx, x_idx]) * tile_dask_array.chunksize
        voxel_max = voxel_min + chunk_shape

        # Map voxel bounding box to physical space
#        chunk_corners = np.array([
#            [voxel_min[2], voxel_min[1], voxel_min[0], 1],
#            [voxel_max[2], voxel_min[1], voxel_min[0], 1],
#            [voxel_min[2], voxel_max[1], voxel_min[0], 1],
#            [voxel_min[2], voxel_min[1], voxel_max[0], 1],
#            [voxel_max[2], voxel_max[1], voxel_min[0], 1],
#            [voxel_min[2], voxel_max[1], voxel_max[0], 1],
#            [voxel_max[2], voxel_min[1], voxel_max[0], 1],
#            [voxel_max[2], voxel_max[1], voxel_max[0], 1],
#        ])
#        print('chunk_corners')
#        print(chunk_corners)
#        physical_corners = np.dot(affine, chunk_corners.T).T[:, :3]

        physical_corners = np.array([
            affine @ np.array([voxel_min[2], voxel_min[1], voxel_min[0], 1]),
            affine @ np.array([voxel_max[2], voxel_min[1], voxel_min[0], 1]),
            affine @ np.array([voxel_min[2], voxel_max[1], voxel_min[0], 1]),
            affine @ np.array([voxel_min[2], voxel_min[1], voxel_max[0], 1]),
            affine @ np.array([voxel_max[2], voxel_max[1], voxel_min[0], 1]),
            affine @ np.array([voxel_min[2], voxel_max[1], voxel_max[0], 1]),
            affine @ np.array([voxel_max[2], voxel_min[1], voxel_max[0], 1]),
            affine @ np.array([voxel_max[2], voxel_max[1], voxel_max[0], 1]),
        ])



#        print('physical corners')
#        print(physical_corners)
        chunk_bbox_min = physical_corners.min(axis=0)[:3]
        chunk_bbox_max = physical_corners.max(axis=0)[:3]

#        print(block_bbox_min)
#        print(chunk_bbox_min)

        # Compute overlapping region in physical space
        overlap_min = np.maximum(block_bbox_min, chunk_bbox_min)
        overlap_max = np.minimum(block_bbox_max, chunk_bbox_max)
        if np.any(overlap_min >= overlap_max):
#            print('NO OVERLAP')
            continue  # No overlap
#        else:
#            print('YES THERE IS OVERLAP')

        # Resample chunk data to the output block's coordinate space
        resampled_chunk = resample_chunk_to_block(
            tile_chunk, affine, overlap_min, overlap_max, block_bbox_min, block_bbox_max, block_shape
        )

    #    print(resampled_chunk.shape)
        # Fuse resampled chunk into the output block
        mask = resampled_chunk > 0
        fused_block[mask] += resampled_chunk[mask]
        weight[mask] += 1

    
    # Normalize fused data
    print('about to fuse')

    fused_block = np.divide(fused_block, weight, out=np.zeros_like(fused_block), where=weight > 0)

    return fused_block


def resample_chunk_to_block(tile_chunk, affine, overlap_min, overlap_max, block_bbox_min, block_bbox_max, block_shape):
    """
    Resample a tile chunk to the output block's coordinate space.

    Parameters:
    - tile_chunk (np.ndarray): The input chunk data from the tile.
    - affine (np.ndarray): The affine transformation of the tile.
    - overlap_min (np.ndarray): Minimum physical coordinates of the overlapping region.
    - overlap_max (np.ndarray): Maximum physical coordinates of the overlapping region.
    - block_bbox_min (np.ndarray): Minimum physical coordinates of the block.
    - block_bbox_max (np.ndarray): Maximum physical coordinates of the block.
    - block_shape (tuple): Shape of the output block.

    Returns:
    - np.ndarray: Resampled data.
    """
    # Compute the overlapping region's voxel coordinates in the tile
    overlap_voxel_min = np.round((overlap_min - affine[:3, 3]) / np.diag(affine[:3, :3])).astype(int)
    overlap_voxel_max = np.round((overlap_max - affine[:3, 3]) / np.diag(affine[:3, :3])).astype(int)

    # Extract the overlapping region from the tile chunk
    overlap_tile = tile_chunk[
        slice(overlap_voxel_min[2], overlap_voxel_max[2]),
        slice(overlap_voxel_min[1], overlap_voxel_max[1]),
        slice(overlap_voxel_min[0], overlap_voxel_max[0]),
    ]

    # Compute the corresponding coordinates in the block
    overlap_block_min = (overlap_min - block_bbox_min) / (block_bbox_max - block_bbox_min) * np.array(block_shape)
    overlap_block_max = (overlap_max - block_bbox_min) / (block_bbox_max - block_bbox_min) * np.array(block_shape)

    # Resample the data using map_coordinates
    coords = np.meshgrid(
        np.linspace(overlap_block_min[0], overlap_block_max[0], overlap_tile.shape[0]),
        np.linspace(overlap_block_min[1], overlap_block_max[1], overlap_tile.shape[1]),
        np.linspace(overlap_block_min[2], overlap_block_max[2], overlap_tile.shape[2]),
        indexing="ij",
    )
    resampled_chunk = map_coordinates(overlap_tile, coords, order=1, mode="constant", cval=0)

    return resampled_chunk




# Paths to input tiles
ome_zarr_paths = snakemake.input.ome_zarr

# Build R-tree index for all chunks in all tiles
rtree_idx, znimgs = build_rtree_index(ome_zarr_paths)


# Parameters for the fused array
fused_shape = (128, 128, 128)  # Full shape of the fused array
chunk_shape = (32,32,32)  # Chunk size of the fused array
global_bbox_min = np.array([-56, -65, -96])  # Minimum physical coordinate of the fused array
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
