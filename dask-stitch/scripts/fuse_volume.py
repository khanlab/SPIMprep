import numpy as np
import dask.array as da
from zarrnii import ZarrNii
from scipy.ndimage import map_coordinates
import os


def compute_fused_volume_shape(ome_zarr_paths, optimized_translations):
    """
    Compute the shape of the fused volume in physical space.

    Parameters:
    - ome_zarr_paths (list of str): List of paths to OME-Zarr datasets.
    - optimized_translations (np.ndarray): Optimized translations for each tile.

    Returns:
    - bbox_min (np.ndarray): Minimum coordinates of the fused volume.
    - bbox_max (np.ndarray): Maximum coordinates of the fused volume.
    - voxel_size (np.ndarray): Voxel size in physical units.
    """
    bbox_min = np.inf * np.ones(3)
    bbox_max = -np.inf * np.ones(3)

    for path, translation in zip(ome_zarr_paths, optimized_translations):
        znimg = ZarrNii.from_path(path)
        affine = znimg.vox2ras.affine

        #HACK fix:
        affine[:3,3] = -1 * np.flip(affine[:3,3])

        img_shape = znimg.darr.shape[1:]  # Exclude the channel/time dimension

        # Compute corrected bounding box with translation
        corners = [
            np.array([0, 0, 0, 1]),
            np.array([img_shape[2], 0, 0, 1]),
            np.array([0, img_shape[1], 0, 1]),
            np.array([img_shape[2], img_shape[1], 0, 1]),
            np.array([0, 0, img_shape[0], 1]),
            np.array([img_shape[2], 0, img_shape[0], 1]),
            np.array([0, img_shape[1], img_shape[0], 1]),
            np.array([img_shape[2], img_shape[1], img_shape[0], 1]),
        ]
        corners_physical = np.dot(affine, np.array(corners).T).T[:, :3] + translation

        # Update bounding box
        bbox_min = np.minimum(bbox_min, corners_physical.min(axis=0))
        bbox_max = np.maximum(bbox_max, corners_physical.max(axis=0))

    voxel_size = np.array((-affine[0,2],-affine[1,1],-affine[2,0]))
    return bbox_min, bbox_max, voxel_size


def resample_tile_to_chunk(tile, affine, translation, chunk_bbox_min, chunk_bbox_max, chunk_shape):
    """
    Resample a tile to the chunk's bounding box in physical space.

    Parameters:
    - tile (np.ndarray): Input tile.
    - affine (np.ndarray): Affine transformation matrix (4x4).
    - translation (np.ndarray): Optimized translation for the tile.
    - chunk_bbox_min (np.ndarray): Minimum physical coordinates of the chunk.
    - chunk_bbox_max (np.ndarray): Maximum physical coordinates of the chunk.
    - chunk_shape (tuple): Shape of the output chunk.

    Returns:
    - np.ndarray: Resampled tile for the chunk.
    """
    coords = np.meshgrid(
        np.linspace(chunk_bbox_min[0], chunk_bbox_max[0], chunk_shape[0]),
        np.linspace(chunk_bbox_min[1], chunk_bbox_max[1], chunk_shape[1]),
        np.linspace(chunk_bbox_min[2], chunk_bbox_max[2], chunk_shape[2]),
        indexing="ij",
    )
    coords = np.stack(coords, axis=-1)  # Shape: (Z, Y, X, 3)

    # Convert physical coordinates to voxel indices
    inverse_affine = np.linalg.inv(affine)
    voxel_coords = np.dot(coords.reshape(-1, 3), inverse_affine[:3, :3].T) + inverse_affine[:3, 3] - translation
    voxel_coords = voxel_coords.reshape(chunk_shape + (3,))

    # Interpolate tile values at the transformed voxel coordinates
    return map_coordinates(tile, [voxel_coords[..., dim] for dim in range(3)], order=1, mode="constant")


def process_chunk(block_info, ome_zarr_paths, optimized_translations, bbox_min, voxel_size, chunk_shape):
    """
    Process a single chunk of the fused volume.

    Parameters:
    - block_info (dict): Information about the block being processed.
    - ome_zarr_paths (list of str): List of paths to OME-Zarr datasets.
    - optimized_translations (np.ndarray): Optimized translations for each tile.
    - bbox_min (np.ndarray): Minimum physical coordinates of the fused volume.
    - voxel_size (np.ndarray): Voxel size in physical units.
    - chunk_shape (tuple): Shape of the chunk being processed.

    Returns:
    - np.ndarray: Fused chunk.
    """
    # Determine chunk bounding box in physical space
    chunk_start = np.array(block_info[0]["array-location"][0]) * voxel_size + bbox_min
    chunk_end = chunk_start + np.array(chunk_shape) * voxel_size

    chunk = np.zeros(chunk_shape, dtype=np.float32)
    weight = np.zeros(chunk_shape, dtype=np.float32)

    for path, translation in zip(ome_zarr_paths, optimized_translations):
        znimg = ZarrNii.from_path(path)
        tile = znimg.darr.squeeze().compute()
        affine = znimg.vox2ras.affine
        affine[:3,3] = -1 * np.flip(affine[:3,3])

        # Resample tile to chunk
        resampled_tile = resample_tile_to_chunk(tile, affine, translation, chunk_start, chunk_end, chunk_shape)

        # Fuse by summing intensities and weights
        mask = resampled_tile > 0
        chunk[mask] += resampled_tile[mask]
        weight[mask] += 1

    # Avoid division by zero
    fused_chunk = np.divide(chunk, weight, out=np.zeros_like(chunk), where=weight > 0)
    return fused_chunk


def fuse_volume(ome_zarr_paths, optimized_translations, fused_shape, chunk_shape, bbox_min, voxel_size):
    """
    Fuse all tiles into a single volume.

    Parameters:
    - ome_zarr_paths (list of str): List of paths to OME-Zarr datasets.
    - optimized_translations (np.ndarray): Optimized translations for each tile.
    - fused_shape (tuple): Shape of the final fused volume.
    - chunk_shape (tuple): Shape of each chunk.
    - bbox_min (np.ndarray): Minimum coordinates of the fused volume.
    - voxel_size (np.ndarray): Voxel size in physical units.

    Returns:
    - dask.array: Fused volume.
    """
    # Wrap process_chunk for Dask
    def wrapped_process_chunk(block, block_info=None):
        return process_chunk(
            block_info,
            ome_zarr_paths,
            optimized_translations,
            bbox_min,
            voxel_size,
            chunk_shape,
        )

    # Define Dask array for the fused volume
    fused_volume = da.map_blocks(
        wrapped_process_chunk,
        chunks=chunk_shape,
        dtype=np.float32,
        #shape=fused_shape,
    )

    return fused_volume


# Example usage
ome_zarr_paths = snakemake.input.ome_zarr  # List of input OME-Zarr paths
optimized_translations = np.loadtxt(snakemake.input.optimized_translations, dtype=float)  # Optimized translations

# Compute the fused volume shape
bbox_min, bbox_max, voxel_size = compute_fused_volume_shape(ome_zarr_paths, optimized_translations)
assert all(voxel_size > 0), "Voxel size must be greater than zero."

fused_shape = tuple(np.ceil((bbox_max - bbox_min) / voxel_size).astype(int))
chunk_shape = (64, 64, 64)  # Example chunk size

# Fuse the volume
fused_volume = fuse_volume(ome_zarr_paths, optimized_translations, fused_shape, chunk_shape, bbox_min, voxel_size)

print(fused_volume.shape)

# Save the fused volume
znimg = ZarrNii.from_darr(fused_volume)

znimg.to_ome_zarr(snakemake.output.ome_zarr)
znimg.to_nifti(snakemake.output.nifti)





