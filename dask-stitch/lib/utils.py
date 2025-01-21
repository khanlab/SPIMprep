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
        # Compute bounding boxes for each chunk
        chunk_bboxes = compute_chunk_bounding_boxes(dask_array, affine, tile_index)
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
def process_block(block_data, block_info=None, rtree_idx=None, znimgs=None, affine=None):
    """
    Process a block of the fused array by loading and resampling overlapping chunks.

    Parameters:
    - block_data (np.ndarray): Placeholder array for the output block (filled with zeros).
    - block_info (dict): Information about the block being processed (from Dask).
    - rtree_idx (rtree.index.Index): R-tree index of the tile chunks.
    - znimgs (list): List of ZarrNii objects for the tiles.
    - affine: Affine transform for output image.

    Returns:
    - np.ndarray: The processed block of the fused array.
    """
    # Extract block location and shape
    block_shape = block_data.shape  # Shape of the block in voxels
    block_location = block_info[0]["chunk-location"]  # (z, y, x) block indices
    array_location = block_info[0]["array-location"]  # (z, y, x) block indices

    # Convert block location to voxel coordinates in the global space
    block_voxel_min = np.array(block_location) * block_shape  # Start voxel indices in global space
    block_voxel_max = block_voxel_min + block_shape  # End voxel indices in global space

    #use the array location instead
    block_bbox_min_vox = (array_location[0][0],
                        array_location[1][0],
                        array_location[2][0])
    block_bbox_max_vox = (array_location[0][1],
                        array_location[1][1],
                        array_location[2][1])


    block_bbox_min = affine @ (np.array(block_bbox_min_vox))
    block_bbox_max = affine @ (np.array(block_bbox_max_vox))

    print('min vox to world:')
    print(block_bbox_min_vox)
    print(affine)
    print(block_bbox_min)

    # Convert voxel coordinates to world coordinates
#    block_bbox_min = global_bbox_min + block_voxel_min * target_voxel_spacing
#    block_bbox_max = global_bbox_min + block_voxel_max * target_voxel_spacing

    #block_bbox = (*block_bbox_min, *block_bbox_max)
    block_bbox = (*(block_bbox_min-5), *(block_bbox_max+5)) #this is for the rtree 


    print(f'at block {block_location}')
    print(f'array_location {array_location}')
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
        affine = znimg.affine.matrix
        print(f'intersected with {tile_idx}')
    

        # Extract chunk data
        z_idx, y_idx, x_idx = chunk_id
        tile_chunk = tile_dask_array.blocks[z_idx, y_idx, x_idx].compute()

        # Compute the bounding box of the chunk in physical space
        chunk_shape = tile_chunk.shape
        voxel_min = np.array([z_idx, y_idx, x_idx]) * tile_dask_array.chunksize
        voxel_max = voxel_min + chunk_shape   #TODO check this

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
   #     print("output chunk extents:")
   #     print(block_bbox_min)
   #     print(block_bbox_max)
   #     print("intersecting input chunk extents:")
   #     print(chunk_bbox_min)
   #     print(chunk_bbox_max)

        # Compute overlapping region in physical space
        overlap_min = np.maximum(block_bbox_min, chunk_bbox_min)
        overlap_max = np.minimum(block_bbox_max, chunk_bbox_max)

   #     print("overlapping region")
   #     print(overlap_min)
   #     print(overlap_max)

#        if np.any(overlap_min > overlap_max):
 #           print('NO OVERLAP')
#            continue  # No overlap
  #      else:
  #          print('YES THERE IS OVERLAP')

        resampled_chunk = resample_chunk_to_block_alt1(
            tile_chunk, chunk_bbox_min, chunk_bbox_max, block_bbox_min, block_bbox_max, block_shape
        )

        #this is new for alt2:
        """
        block_coords = [
            np.linspace(block_bbox_min[i], block_bbox_max[i], block_shape[i])
            for i in range(3)
        ]
        output_coords = np.meshgrid(block_coords[2], block_coords[1], block_coords[0], indexing="ij")
        resampled_chunk = resample_chunk_to_block_alt2(tile_chunk, affine, output_coords)
        """
        #--

        # Resample chunk data to the output block's coordinate space
       # resampled_chunk = resample_chunk_to_block(
       #     tile_chunk, affine, overlap_min, overlap_max, block_bbox_min, block_bbox_max, block_shape
       # )

    #    print(resampled_chunk.shape)
        # Fuse resampled chunk into the output block
        mask = resampled_chunk > 0
        fused_block[mask] += resampled_chunk[mask]
        weight[mask] += 1

    
    # Normalize fused data
    print('about to fuse')

    fused_block = np.divide(fused_block, weight, out=np.zeros_like(fused_block), where=weight > 0)

    return fused_block


def resample_chunk_to_block_alt2(tile_chunk, affine, output_coords):
    """
    Resample input chunk data to the output block's coordinate space using an affine transformation.

    Parameters:
    - tile_chunk (np.ndarray): Input chunk data (3D array).
    - affine (np.ndarray): 4x4 affine matrix (voxel-to-world coordinates).
    - output_coords (tuple): Meshgrid of physical coordinates for the output block:
        (output_x, output_y, output_z).

    Returns:
    - np.ndarray: Resampled data for the output block.
    """
    # Compute the inverse affine to map world to voxel coordinates
    affine_inv = np.linalg.inv(affine)

    # Map output physical coordinates to input chunk voxel coordinates
    input_coords = [
        affine_inv[i, :3] @ np.stack(output_coords, axis=-1).reshape(-1, 3).T + affine_inv[i, 3]
        for i in range(3)
    ]

    # Reshape input_coords back to the shape of output_coords
    input_coords = [coord.reshape(output_coords[0].shape) for coord in input_coords]

    # Interpolate the chunk data onto the output block
    resampled_chunk = map_coordinates(tile_chunk, input_coords, order=2, mode="constant", cval=0)

    return resampled_chunk


def resample_chunk_to_block_alt1(tile_chunk, chunk_bbox_min, chunk_bbox_max, block_bbox_min, block_bbox_max, block_shape):
    """
    Resample tile_chunk (input chunk) data to the output block's coordinate space.

    Parameters:
    - tile_chunk (np.ndarray): Input chunk data (3D array).
    - chunk_bbox_min (np.ndarray): Physical min coordinates of the input chunk.
    - chunk_bbox_max (np.ndarray): Physical max coordinates of the input chunk.
    - block_bbox_min (np.ndarray): Physical min coordinates of the output block.
    - block_bbox_max (np.ndarray): Physical max coordinates of the output block.
    - block_shape (tuple): Shape of the output block.

    Returns:
    - np.ndarray: Resampled data for the output block.
    """
    # Create output coordinates for the block (linspace in physical space)
    block_coords = [
        np.linspace(block_bbox_min[i], block_bbox_max[i], block_shape[i])
        for i in range(3)
    ]

    # Meshgrid for output block coordinates (physical space)
    output_coords = np.meshgrid(
        block_coords[2],  # x-coordinates
        block_coords[1],  # y-coordinates
        block_coords[0],  # z-coordinates
        indexing="ij",
    )

    # Map output block coordinates to input chunk voxel coordinates
    input_coords = [
        (output_coords[i] - chunk_bbox_min[i]) /
        (chunk_bbox_max[i] - chunk_bbox_min[i]) * (tile_chunk.shape[i] - 1)
        for i in range(3)
    ]

    # Interpolate tile_chunk data onto the output block's coordinate space
    resampled_chunk = map_coordinates(tile_chunk, input_coords, order=1, mode="constant", cval=0)

    return resampled_chunk



