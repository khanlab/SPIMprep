import numpy as np
from zarrnii import ZarrNii
from scipy.fft import fftn, ifftn
from scipy.ndimage import map_coordinates


def phase_correlation(img1, img2):
    """
    Compute the phase correlation between two 3D images to find the translation offset.

    Parameters:
    - img1 (np.ndarray): First image (3D array).
    - img2 (np.ndarray): Second image (3D array).

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


def compute_pairwise_correlation(ome_zarr_paths, overlapping_pairs, output_shape=(64, 64, 64)):
    """
    Compute the optimal offset for each pair of overlapping tiles.

    Parameters:
    - ome_zarr_paths (list of str): List of paths to OME-Zarr datasets.
    - overlapping_pairs (list of tuples): List of overlapping tile indices.
    - output_shape (tuple): Shape of the resampled bounding box.

    Returns:
    - np.ndarray: Array of offsets for each pair (N, 3) where N is the number of pairs.
    """
    offsets = []

    for i, j in overlapping_pairs:
        # Load the two images and their affines
        znimg1 = ZarrNii.from_path(ome_zarr_paths[i])
        znimg2 = ZarrNii.from_path(ome_zarr_paths[j])

        img1 = znimg1.darr.squeeze().compute()
        img2 = znimg2.darr.squeeze().compute()

        affine1 = znimg1.vox2ras.affine
        affine2 = znimg2.vox2ras.affine

        #HACK FIX
        affine1[:3,3] = -1 * np.flip(affine1[:3,3])
        affine2[:3,3] = -1 * np.flip(affine2[:3,3])

        # Compute the overlapping bounding box in physical space
        bbox1_min = affine1[:3, 3]
        bbox1_max = bbox1_min + np.dot(affine1[:3, :3], img1.shape[::-1])
        bbox2_min = affine2[:3, 3]
        bbox2_max = bbox2_min + np.dot(affine2[:3, :3], img2.shape[::-1])

        bbox_min = np.maximum(bbox1_min, bbox2_min)
        bbox_max = np.minimum(bbox1_max, bbox2_max)

        # Resample images to the overlapping bounding box
        resampled_img1 = resample_to_bounding_box(img1, affine1, bbox_min, bbox_max, output_shape)
        resampled_img2 = resample_to_bounding_box(img2, affine2, bbox_min, bbox_max, output_shape)

        # Compute phase correlation on the resampled overlapping region
        offset = phase_correlation(resampled_img1, resampled_img2)
        offsets.append(offset)

    return np.array(offsets)


# Example usage
overlapping_pairs = np.loadtxt(snakemake.input.pairs, dtype=int).tolist()  # Overlapping pairs
ome_zarr_paths = snakemake.input.ome_zarr  # List of OME-Zarr paths

# Compute pairwise offsets
offsets = compute_pairwise_correlation(ome_zarr_paths, overlapping_pairs)

# Save results
np.savetxt(snakemake.output.offsets, offsets, fmt="%.6f")


