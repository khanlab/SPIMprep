import numpy as np
from templateflow import api as tflow
import dask.array as da
import math
from zarrnii import ZarrNii

def create_test_dataset(template="MNI152NLin2009cAsym", res=2, grid_shape=(3, 4), overlap=8, random_seed=42, final_chunks=(32, 32, 1)):
    """
    Create a low-resolution test dataset for tile-based stitching.

    Parameters:
    - template (str): TemplateFlow template name (default: MNI152NLin2009cAsym).
    - res (int): Desired resolution in mm (default: 2mm).
    - grid_shape (tuple): Shape of the tiling grid (e.g., (3, 4) for 3x4 grid in X-Y).
    - overlap (int): Overlap between tiles in X-Y plane in voxels (default: 8).
    - random_seed (int): Seed for reproducible random offsets.
    - final_chunks (tuple): Desired chunks for final tiles.

    Returns:
    - znimg (ZarrNii): ZarrNii object containing the tiles.
    - translations (np.ndarray): Array of random offsets for each tile.
    """
    import math

    # Seed the random number generator
    np.random.seed(random_seed)

    # Download template and load as a ZarrNii object
    template_path = tflow.get(template, resolution=res, desc=None, suffix="T1w")
    znimg = ZarrNii.from_path(template_path)
    img_data = znimg.darr.squeeze()

    # Original image shape
    img_shape = np.array(img_data.shape)  # (Z, Y, X)

    # Keep Z dimension intact, calculate X and Y tile sizes
    z_dim, y_dim, x_dim = img_shape
    x_tile_size = math.ceil((x_dim + overlap * (grid_shape[1] - 1)) / grid_shape[1])
    y_tile_size = math.ceil((y_dim + overlap * (grid_shape[0] - 1)) / grid_shape[0])
    tile_shape = (z_dim, y_tile_size, x_tile_size)

    # Calculate required padding to ensure X and Y dimensions are divisible by grid shape
    padded_x = x_tile_size * grid_shape[1] - overlap * (grid_shape[1] - 1)
    padded_y = y_tile_size * grid_shape[0] - overlap * (grid_shape[0] - 1)

    padding = (
        (0, 0),  # No padding in Z
        (0, int(max(padded_y - y_dim, 0))),
        (0, int(max(padded_x - x_dim, 0))),
    )

    print('padding')
    print(padding)
    print(img_data.shape)

    # Pad image if needed
    if any(p[1] > 0 for p in padding):
        img_data = da.pad(img_data, padding, mode="constant", constant_values=0)

    # Create tiles
    tiles = []
    translations = []
    for y in range(grid_shape[0]):
        for x in range(grid_shape[1]):
            # Calculate tile start indices
            y_start = y * (y_tile_size - overlap)
            x_start = x * (x_tile_size - overlap)

            # Extract tile
            tile = img_data[:, y_start:y_start + y_tile_size, x_start:x_start + x_tile_size]
            tiles.append(tile)

            # Add random offset -- NOT ACTUALLY BEING APPLIED TO SAMPLING HERE!
            offset = np.random.uniform(-5, 5, size=3)  # Random 3D offsets 
            translations.append((0, y_start, x_start) + offset)

    # Convert tiles to a Dask array
    tiles = da.stack([tile.rechunk(chunks=final_chunks) for tile in tiles])

    # Save back into ZarrNii object
    znimg.darr = tiles
    translations = np.array(translations)

    return znimg, translations




test_znimg, test_translations = create_test_dataset(grid_shape=snakemake.params.grid_shape)
print(test_znimg.darr.shape)
test_znimg.to_ome_zarr(snakemake.output.ome_zarr)
np.save(snakemake.output.translations_npy,test_translations)


