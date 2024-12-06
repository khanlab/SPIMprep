import numpy as np

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
        znimg = ZarrNii.from_path(path)
        affine = znimg.vox2ras.affine  # 4x4 matrix
        tile_shape = znimg.darr.shape[1:]

        # Compute physical bounding box using affine
        corners = [
            np.array([0, 0, 0, 1]),
            np.array([tile_shape[2], 0, 0, 1]),
            np.array([0, tile_shape[1], 0, 1]),
            np.array([tile_shape[2], tile_shape[1], 0, 1]),
            np.array([0, 0, tile_shape[0], 1]),
            np.array([tile_shape[2], 0, tile_shape[0], 1]),
            np.array([0, tile_shape[1], tile_shape[0], 1]),
            np.array([tile_shape[2], tile_shape[1], tile_shape[0], 1]),
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


overlapping_pairs = find_overlapping_pairs(snakemake.input)
np.savetxt(snakemake.output.txt,np.array(overlapping_pairs),fmt='%d')

