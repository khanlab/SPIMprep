import numpy as np
from zarrnii import ZarrNii
from scipy.fft import fftn, ifftn
from scipy.ndimage import center_of_mass


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


def compute_pairwise_correlation(ome_zarr_paths, overlapping_pairs):
    """
    Compute the optimal offset for each pair of overlapping tiles.

    Parameters:
    - ome_zarr_paths (list of str): List of paths to OME-Zarr datasets.
    - overlapping_pairs (list of tuples): List of overlapping tile indices.

    Returns:
    - np.ndarray: Array of offsets for each pair (N, 3) where N is the number of pairs.
    """
    offsets = []

    for i, j in overlapping_pairs:
        # Load the two images
        znimg1 = ZarrNii.from_path(ome_zarr_paths[i])
        znimg2 = ZarrNii.from_path(ome_zarr_paths[j])

        img1 = znimg1.darr.squeeze().compute()
        img2 = znimg2.darr.squeeze().compute()

        # Compute phase correlation
        offset = phase_correlation(img1, img2)
        offsets.append(offset)

    return np.array(offsets)


# Example usage
overlapping_pairs = np.loadtxt(snakemake.input.pairs, dtype=int).tolist()  # Overlapping pairs
ome_zarr_paths = snakemake.input.ome_zarr  # List of OME-Zarr paths

# Compute pairwise offsets
offsets = compute_pairwise_correlation(ome_zarr_paths, overlapping_pairs)

# Save results
np.savetxt(snakemake.output.offsets, offsets, fmt="%.6f")

