import numpy as np
from scipy.optimize import least_squares
from zarrnii import ZarrNii


def global_optimization(ome_zarr_paths, overlapping_pairs, pairwise_offsets):
    """
    Perform global optimization to adjust translations for all tiles.

    Parameters:
    - ome_zarr_paths (list of str): List of paths to OME-Zarr datasets.
    - overlapping_pairs (list of tuples): List of overlapping tile indices ((i, j)).
    - pairwise_offsets (np.ndarray): Array of pairwise offsets (N, 3), where N is the number of pairs.

    Returns:
    - np.ndarray: Optimized global translations of shape (T, 3).
    """
    # Number of tiles is the number of OME-Zarr paths
    num_tiles = len(ome_zarr_paths)

    # Initial translations (start with identity translation: no offsets)
    initial_translations = np.zeros((num_tiles, 3))

    # Flatten initial translations for optimization
    x0 = initial_translations.flatten()

    def objective(x):
        """
        Compute the residuals for global optimization.

        Parameters:
        - x (np.ndarray): Flattened translations array (T * 3,).

        Returns:
        - np.ndarray: Residuals for least-squares optimization.
        """
        translations = x.reshape((num_tiles, 3))
        residuals = []

        for (i, j), offset in zip(overlapping_pairs, pairwise_offsets):
            # Residual is the difference between the predicted and actual offset
            predicted_offset = translations[j] - translations[i]
            residuals.append(predicted_offset - offset)

        return np.concatenate(residuals)

    # Perform least-squares optimization
    result = least_squares(objective, x0)

    # Reshape result back to (T, 3)
    optimized_translations = result.x.reshape((num_tiles, 3))

    return optimized_translations


# Example usage
overlapping_pairs = np.loadtxt(snakemake.input.pairs, dtype=int).tolist()  # Overlapping pairs
pairwise_offsets = np.loadtxt(snakemake.input.offsets, dtype=float)  # Pairwise offsets
ome_zarr_paths = snakemake.input.ome_zarr  # List of OME-Zarr paths

# Perform global optimization
optimized_translations = global_optimization(ome_zarr_paths, overlapping_pairs, pairwise_offsets)

# Save results
np.savetxt(snakemake.output.optimized_translations, optimized_translations, fmt="%.6f")

