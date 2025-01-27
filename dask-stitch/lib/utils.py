import numpy as np
import dask.array as da
from zarrnii import ZarrNii, AffineTransform
from scipy.ndimage import map_coordinates
from rtree import index

#import dask
#dask.config.set(scheduler='synchronous')  # overwrite default with single-threaded scheduler


def find_overlapping_pairs(ome_zarr_paths):
    """
    Identify overlapping tile pairs based on their physical offsets.

    Parameters:
    - ome_zarr_paths (list of str): List of paths to OME-Zarr datasets.

    Returns:
    - List of tuples: Each tuple is a pair of overlapping tile indices ((i, j)).
    """
    from zarrnii import ZarrNii

    # Read physical transformations and calculate bounding boxes
    bounding_boxes = []
    for path in ome_zarr_paths:
        znimg = ZarrNii.from_ome_zarr(path)

        tile_shape = znimg.darr.shape[1:]
 
        affine = np.eye(4,4)
        affine[:3, :3] = np.diag(znimg.get_zooms(axes_order=znimg.axes_order))  # Set the zooms (scaling factors) along the diagonal
        affine[:3, 3] = znimg.get_origin(axes_order=znimg.axes_order)  # Set the translation (origin)
 
        # Compute physical bounding box using affine
        corners = [
            np.array([0, 0, 0, 1]),
            np.array([tile_shape[0], 0, 0, 1]),
            np.array([0, tile_shape[1], 0, 1]),
            np.array([tile_shape[0], tile_shape[1], 0, 1]),
            np.array([0, 0, tile_shape[0], 1]),
            np.array([tile_shape[0], 0, tile_shape[2], 1]),
            np.array([0, tile_shape[1], tile_shape[2], 1]),
            np.array([tile_shape[0], tile_shape[1], tile_shape[2], 1]),
        ]
        corners_physical = np.dot(affine, np.array(corners).T).T[:, :3]  # Drop homogeneous coordinate
        bbox_min = corners_physical.min(axis=0)
        bbox_max = corners_physical.max(axis=0)

        bounding_boxes.append((bbox_min, bbox_max))

    # Find overlapping pairs
    overlapping_pairs = []
    for i, (bbox1_min, bbox1_max) in enumerate(bounding_boxes):
        for j, (bbox2_min, bbox2_max) in enumerate(bounding_boxes):
            if i >= j:
                continue  # Avoid duplicate pairs and self-comparison

            # Check for overlap in all dimensions
            overlap = all(
                bbox1_min[d] < bbox2_max[d] and bbox1_max[d] > bbox2_min[d]
                for d in range(3)
            )
            if overlap:
                overlapping_pairs.append((i, j))

    return overlapping_pairs


def compute_chunk_bounding_boxes(dask_array, zooms, origin, tile_index):
    """
    Compute the bounding boxes of each chunk in a Dask array in physical space.

    Parameters:
    - dask_array (dask.array): The input tile as a Dask array.
    - zooms (np.ndarray): Voxel spacing 
    - origin (np.ndarray): Offset or translation of origin.
    - tile_index (int): Index of the tile for tracking.

    Returns:
    - List of bounding boxes in physical space.
      Each bounding box is (chunk_id, (min_x, min_y, min_z, max_x, max_y, max_z)).
    """

    chunk_bboxes = []
    chunk_shapes = dask_array.chunks  # Tuple of chunk sizes along each dimension

    #make affine from zooms and origin
    affine = np.eye(4, dtype=float)
    affine[:3, :3] = np.diag(zooms)  # Set the zooms (scaling factors) along the diagonal
    affine[:3, 3] = origin  # Set the translation (origin)

    # Iterate through all chunk indices
    for x_idx, x_start in enumerate(np.cumsum((0,) + chunk_shapes[0][:-1])):
        x_end = x_start + chunk_shapes[0][x_idx]

        for y_idx, y_start in enumerate(np.cumsum((0,) + chunk_shapes[1][:-1])):
            y_end = y_start + chunk_shapes[1][y_idx]

            for z_idx, z_start in enumerate(np.cumsum((0,) + chunk_shapes[2][:-1])):
                z_end = z_start + chunk_shapes[2][z_idx]

                # Voxel coordinates for the chunk
                voxel_min = np.array([x_start, y_start, z_start])
                voxel_max = np.array([x_end, y_end, z_end])

                # Convert voxel coordinates to physical space using affine
                physical_corners = np.array([
                    affine @ np.array([voxel_min[0], voxel_min[1], voxel_min[2], 1]),
                    affine @ np.array([voxel_max[0], voxel_min[1], voxel_min[2], 1]),
                    affine @ np.array([voxel_min[0], voxel_max[1], voxel_min[2], 1]),
                    affine @ np.array([voxel_min[0], voxel_min[1], voxel_max[2], 1]),
                    affine @ np.array([voxel_max[0], voxel_max[1], voxel_min[2], 1]),
                    affine @ np.array([voxel_min[0], voxel_max[1], voxel_max[2], 1]),
                    affine @ np.array([voxel_max[0], voxel_min[1], voxel_max[2], 1]),
                    affine @ np.array([voxel_max[0], voxel_max[1], voxel_max[2], 1]),
                ])

                physical_bbox_min = physical_corners.min(axis=0)[:3]
                physical_bbox_max = physical_corners.max(axis=0)[:3]

                # Store the bounding box with a unique chunk_id
                chunk_id = (x_idx, y_idx, z_idx)
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

        #TODO: get_origin and get_zooms should specify an axes_order..

        # Compute bounding boxes for each chunk
        chunk_bboxes = compute_chunk_bounding_boxes(dask_array, zooms=znimg.get_zooms(axes_order=znimg.axes_order),origin=znimg.get_origin(axes_order=znimg.axes_order), tile_index=tile_index)
        # Insert bounding boxes into the R-tree
        for chunk_id, (physical_bbox_min, physical_bbox_max, tile_idx) in chunk_bboxes:
            #print(f'{chunk_id}, {physical_bbox_min} to {physical_bbox_max} in tile {tile_idx}')
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



def resample_chunk_to_block(tile_chunk, chunk_bbox_min, chunk_bbox_max, block_bbox_min, block_bbox_max, block_shape):
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
        block_coords[0],  # x-coordinates
        block_coords[1],  # y-coordinates
        block_coords[2],  # z-coordinates
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




# Define the process_block function
def process_block(block_data, block_info=None, rtree_idx=None, znimgs=None, ref_znimg=None):
    """
    Process a block of the fused array by loading and resampling overlapping chunks.

    Parameters:
    - block_data (np.ndarray): Placeholder array for the output block (filled with zeros).
    - block_info (dict): Information about the block being processed (from Dask).
    - rtree_idx (rtree.index.Index): R-tree index of the tile chunks.
    - znimgs (list): List of ZarrNii objects for the tiles.
    - ref_znimg (zarrnii.ZarrNii): ZarrNii object for reference (output) blocks.

    Returns:
    - np.ndarray: The processed block of the fused array.
    """
    # Extract block location and shape
    block_shape = block_data.shape  # Shape of the block in voxels
    block_location = block_info[0]["chunk-location"]  # (x, y, z) block indices
    array_location = block_info[0]["array-location"]  # (x, y, z) block indices

    #use the array location instead
    block_bbox_min_vox = (array_location[0][0],
                        array_location[1][0],
                        array_location[2][0])
    block_bbox_max_vox = (array_location[0][1],
                        array_location[1][1],
                        array_location[2][1])

    #make affine from zooms and origin - so that any axis reordering is ignored
    affine = AffineTransform.identity()
    affine[:3, :3] = np.diag(ref_znimg.get_zooms(axes_order=ref_znimg.axes_order))  # Set the zooms (scaling factors) along the diagonal
    affine[:3, 3] = ref_znimg.get_origin(axes_order=ref_znimg.axes_order)  # Set the translation (origin)
    
    

    block_bbox_min = affine @ (np.array(block_bbox_min_vox))
    block_bbox_max = affine @ (np.array(block_bbox_max_vox))

#    print('min vox to world:')
#    print(block_bbox_min_vox)
#    print(affine)
#    print(block_bbox_min)


    #block_bbox = (*block_bbox_min, *block_bbox_max)
    block_bbox = (*(block_bbox_min), *(block_bbox_max)) #this is for the rtree 


#    print(f'at block {block_location}')
#    print(f'array_location {array_location}')
#    print(f'bbox: {block_bbox}')


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
#        print(f'intersected with {tile_idx}')
 
        affine = AffineTransform.identity()
        affine[:3, :3] = np.diag(znimg.get_zooms(axes_order=znimg.axes_order))  # Set the zooms (scaling factors) along the diagonal
        affine[:3, 3] = znimg.get_origin(axes_order=znimg.axes_order)  # Set the translation (origin)
    

        # Extract chunk data
        x_idx, y_idx, z_idx = chunk_id
        tile_chunk = tile_dask_array.blocks[x_idx, y_idx, z_idx].compute()

        # Compute the bounding box of the chunk in physical space
        chunk_shape = tile_chunk.shape
        voxel_min = np.array([x_idx, y_idx, z_idx]) * tile_dask_array.chunksize
        voxel_max = voxel_min + chunk_shape   #TODO check this

        physical_corners = np.array([
                    affine @ np.array([voxel_min[0], voxel_min[1], voxel_min[2], 1]),
                    affine @ np.array([voxel_max[0], voxel_min[1], voxel_min[2], 1]),
                    affine @ np.array([voxel_min[0], voxel_max[1], voxel_min[2], 1]),
                    affine @ np.array([voxel_min[0], voxel_min[1], voxel_max[2], 1]),
                    affine @ np.array([voxel_max[0], voxel_max[1], voxel_min[2], 1]),
                    affine @ np.array([voxel_min[0], voxel_max[1], voxel_max[2], 1]),
                    affine @ np.array([voxel_max[0], voxel_min[1], voxel_max[2], 1]),
                    affine @ np.array([voxel_max[0], voxel_max[1], voxel_max[2], 1]),
                ])


        chunk_bbox_min = physical_corners.min(axis=0)[:3]
        chunk_bbox_max = physical_corners.max(axis=0)[:3]

        resampled_chunk = resample_chunk_to_block(
            tile_chunk, chunk_bbox_min, chunk_bbox_max, block_bbox_min, block_bbox_max, block_shape
        )


        # Fuse resampled chunk into the output block
        mask = resampled_chunk > 0
        fused_block[mask] += resampled_chunk[mask]
        weight[mask] += 1


#    print(f"fusing output block {block_info[0]['chunk-location']}")


    #TODO: commenting this for now just to evaluate bbox intersections
    fused_block = np.divide(fused_block, weight, out=np.zeros_like(fused_block), where=weight > 0)

    return fused_block


def fuse_output_blocks(znimg,rtree_idx, tile_znimgs):

    output_shape = znimg.darr.shape

    # Apply map_blocks to process the fused array --- could consider putting this into zarrnii as a member function.. 
    processed_dask_array = znimg.darr.map_blocks(
        process_block,
        dtype=np.float32,
        rtree_idx=rtree_idx,
        znimgs=tile_znimgs,
        ref_znimg=znimg
    )

    znimg.darr = processed_dask_array.reshape(1,output_shape[0],output_shape[1],output_shape[2])

    return znimg

 

