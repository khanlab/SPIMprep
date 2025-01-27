import numpy as np
from scipy.optimize import least_squares

from scipy.sparse import csr_matrix
from scipy.sparse.linalg import lsqr

def calculate_global_translations(pairs, pairwise_offsets, num_tiles):
    """
    Calculate global translations from pairwise offsets.

    Parameters:
    - pairs (np.ndarray): Array of shape (N, 2) containing indices of overlapping tile pairs.
    - pairwise_offsets (np.ndarray): Array of shape (N, 3) containing pairwise offsets for each pair.
    - num_tiles (int): Total number of tiles.

    Returns:
    - np.ndarray: Global translations for each tile of shape (num_tiles, 3).
    """
    # Number of pairwise offsets
    num_pairs = pairs.shape[0]

    # Initialize the sparse matrix A and vector b
    data = []
    row_indices = []
    col_indices = []
    b = np.zeros((num_pairs * 3,))

    for k, (i, j) in enumerate(pairs):
        # Each pair contributes to three equations (x, y, z components)
        for d in range(3):  # x=0, y=1, z=2
            row = 3 * k + d

            # T_j - T_i = O_ij
            data.append(-1)
            row_indices.append(row)
            col_indices.append(3 * i + d)  # T_i[d]

            data.append(1)
            row_indices.append(row)
            col_indices.append(3 * j + d)  # T_j[d]

            # Right-hand side
            b[row] = pairwise_offsets[k, d]

    # Convert to sparse matrix
    A = csr_matrix((data, (row_indices, col_indices)), shape=(num_pairs * 3, num_tiles * 3))

    # Anchor the first tile (T_0 = [0, 0, 0])
    anchor_rows = np.zeros((3, num_tiles * 3))
    for d in range(3):
        anchor_rows[d, d] = 1
    A = csr_matrix(np.vstack([A.toarray(), anchor_rows]))
    b = np.hstack([b, [0, 0, 0]])

    # Solve the linear system using least squares
    x = lsqr(A, b)[0]

    # Reshape the result into (num_tiles, 3)
    global_translations = x.reshape((num_tiles, 3))

    return global_translations


# Example usage
overlapping_pairs = np.loadtxt(snakemake.input.pairs, dtype=int)  # Overlapping pairs
pairwise_offsets = np.loadtxt(snakemake.input.offsets, dtype=float)  # Pairwise offsets
n_tiles = snakemake.params.n_tiles

# Perform global optimization
optimized_translations = calculate_global_translations(overlapping_pairs, pairwise_offsets, n_tiles)

# Save results
np.savetxt(snakemake.output.optimized_translations, optimized_translations, fmt="%.6f")

