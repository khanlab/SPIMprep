import numpy as np
from zarrnii import ZarrNii
from scipy.fft import fftn, ifftn
from scipy.ndimage import map_coordinates
import os


def phase_correlation(img1, img2, diagnostics_dir=None, pair_index=None):
    """
    Compute the phase correlation between two 3D images to find the translation offset.

    Parameters:
    - img1 (np.ndarray): First image (3D array).
    - img2 (np.ndarray): Second image (3D array).
    - diagnostics_dir (str): Directory to save diagnostic outputs (optional).
    - pair_index (int): Index of the pair being processed (for naming diagnostics).

    Returns:
    - np.ndarray: Offset vector [z_offset, y_offset, x_offset].
    """
    # Compute the Fourier transforms
    fft1 = fftn(img1)
    fft2 = fftn(img2)

    # Compute the cross-power spectrum
    cross_power = fft1 * np.conj(fft2)
    cross_power /= np.abs(cross_power)  # Normalize

    # Inverse Fourier transform to get correlation map
    correlation = np.abs(ifftn(cross_power))

    # Save the correlation map for diagnostics
    if diagnostics_dir and pair_index is not None:
        np.save(os.path.join(diagnostics_dir, f"correlation_map_pair_{pair_index}.npy"), correlation)

    # Find the peak in the correlation map
    peak = np.unravel_index(np.argmax(correlation), correlation.shape)

    # Convert peak index to an offset
    shifts = np.array(peak, dtype=float)
    for dim, size in enumerate(correlation.shape):
        if shifts[dim] > size // 2:
            shifts[dim] -= size

    return shifts


def resample_to_bounding_box(img, affine, bbox_min, bbox_max, output_shape):
    """
    Resample the image to the specified bounding box in physical space.

    Parameters:
    - img (np.ndarray): Input image.
    - affine (np.ndarray): Affine transformation matrix (4x4).
    - bbox_min (np.ndarray): Minimum physical coordinates of the bounding box.
    - bbox_max (np.ndarray): Maximum physical coordinates of the bounding box.
    - output_shape (tuple): Desired shape of the output image.

    Returns:
    - np.ndarray: Resampled image.
    """
    coords = np.meshgrid(
        np.linspace(bbox_min[0], bbox_max[0], output_shape[0]),
        np.linspace(bbox_min[1], bbox_max[1], output_shape[1]),
        np.linspace(bbox_min[2], bbox_max[2], output_shape[2]),
        indexing="ij",
    )
    coords = np.stack(coords, axis=-1)  # Shape: (Z, Y, X, 3)

    # Convert physical coordinates to voxel indices
    inverse_affine = np.linalg.inv(affine)
    voxel_coords = np.dot(coords.reshape(-1, 3), inverse_affine[:3, :3].T) + inverse_affine[:3, 3]
    voxel_coords = voxel_coords.reshape(output_shape + (3,))

    # Interpolate image values at the transformed voxel coordinates
    return map_coordinates(img, [voxel_coords[..., dim] for dim in range(3)], order=1, mode="constant")


def compute_corrected_bounding_box(affine, img_shape):
    """
    Compute the corrected bounding box in physical space, accounting for negative physical dimensions.

    Parameters:
    - affine (np.ndarray): 4x4 affine matrix.
    - img_shape (tuple): Shape of the image (Z, Y, X).

    Returns:
    - bbox_min (np.ndarray): Minimum physical coordinates.
    - bbox_max (np.ndarray): Maximum physical coordinates.
    """
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

    # Transform voxel corners to physical space
    corners_physical = np.dot(affine, np.array(corners).T).T[:, :3]

    # Find corrected min and max
    bbox_min = np.min(corners_physical, axis=0)
    bbox_max = np.max(corners_physical, axis=0)

    return bbox_min, bbox_max


def compute_pairwise_correlation(ome_zarr_paths, overlapping_pairs, output_shape=(64, 64, 64), diagnostics_dir=None):
    """
    Compute the optimal offset for each pair of overlapping tiles.

    Parameters:
    - ome_zarr_paths (list of str): List of paths to OME-Zarr datasets.
    - overlapping_pairs (list of tuples): List of overlapping tile indices.
    - output_shape (tuple): Shape of the resampled bounding box.
    - diagnostics_dir (str): Directory to save diagnostic outputs (optional).

    Returns:
    - np.ndarray: Array of offsets for each pair (N, 3) where N is the number of pairs.
    """
    offsets = []

    if diagnostics_dir:
        os.makedirs(diagnostics_dir, exist_ok=True)

    for pair_index, (i, j) in enumerate(overlapping_pairs):
        # Load the two images and their affines
        znimg1 = ZarrNii.from_path(ome_zarr_paths[i])
        znimg2 = ZarrNii.from_path(ome_zarr_paths[j])

        img1 = znimg1.darr.squeeze().compute()
        img2 = znimg2.darr.squeeze().compute()

        affine1 = znimg1.vox2ras.affine
        affine2 = znimg2.vox2ras.affine
        #HACK FIX
#        affine1[:3,3] = -1 * np.flip(affine1[:3,3])
#        affine2[:3,3] = -1 * np.flip(affine2[:3,3])



        # Compute the corrected bounding boxes
        bbox1_min, bbox1_max = compute_corrected_bounding_box(affine1, img1.shape)
        bbox2_min, bbox2_max = compute_corrected_bounding_box(affine2, img2.shape)

        

        print(f'bbox1, from tile: {i}')
        print(bbox1_min)
        print(bbox1_max)
        print(f'bbox2, from tile: {j}')
        print(bbox2_min)
        print(bbox2_max)

        bbox_min = np.maximum(bbox1_min, bbox2_min)
        bbox_max = np.minimum(bbox1_max, bbox2_max)

        print('overlapping bbox')
        print(bbox_min)
        print(bbox_max)
        # Save bounding box for diagnostics
        if diagnostics_dir:
            np.save(os.path.join(diagnostics_dir, f"bounding_box_pair_{pair_index}.npy"), np.array([bbox_min, bbox_max]))

        # Resample images to the overlapping bounding box
        resampled_img1 = resample_to_bounding_box(img1, affine1, bbox_min, bbox_max, output_shape)
        resampled_img2 = resample_to_bounding_box(img2, affine2, bbox_min, bbox_max, output_shape)

        # Save resampled images for diagnostics
        if diagnostics_dir:
            np.save(os.path.join(diagnostics_dir, f"resampled_img1_pair_{pair_index}.npy"), resampled_img1)
            np.save(os.path.join(diagnostics_dir, f"resampled_img2_pair_{pair_index}.npy"), resampled_img2)

        # Compute phase correlation on the resampled overlapping region
        offset = phase_correlation(resampled_img1, resampled_img2, diagnostics_dir=diagnostics_dir, pair_index=pair_index)
        offsets.append(offset)

    return np.array(offsets)


# Example usage
overlapping_pairs = np.loadtxt(snakemake.input.pairs, dtype=int).tolist()  # Overlapping pairs
ome_zarr_paths = snakemake.input.ome_zarr  # List of OME-Zarr paths
diagnostics_dir = snakemake.output.diagnostics_dir  # Directory for diagnostics

# Compute pairwise offsets
offsets = compute_pairwise_correlation(ome_zarr_paths, overlapping_pairs, diagnostics_dir=diagnostics_dir)

# Save results
np.savetxt(snakemake.output.offsets, offsets, fmt="%.6f")


