import numpy as np
from templateflow import api as tflow
import dask.array as da
from zarrnii import ZarrNii


out_ome_zarr=snakemake.output.ome_zarr
out_translations_npy=snakemake.output.translations_npy

def create_test_dataset(template="MNI152NLin2009cAsym", res=2, tile_shape=(32,32, 32), final_chunks=(1,32,32,1),overlap=8, random_seed=42):
    """
    Create a low-resolution test dataset for tile-based stitching.

    Parameters:
    - template (str): TemplateFlow template name (default: MNI152NLin2009cAsym).
    - res (int): Desired resolution in mm (default: 2mm).
    - tile_shape (tuple): Shape of each tile in voxels (default: (64, 64, 64)).
    - overlap (int): Overlap between tiles in voxels (default: 16).
    - random_seed (int): Seed for reproducible random offsets.

    Returns:
    - tiles (dask.array): TxZxYxX Dask array of tiles.
    - translations (np.ndarray): Array of random offsets for each tile.
    """
    # Seed the random number generator
    np.random.seed(random_seed)

    # Download template and load as a Numpy array
    template_path = tflow.get(template, resolution=res, desc=None,suffix="T1w")
    znimg = ZarrNii.from_path(template_path)
    print(znimg)
    print(znimg.darr)
    img_data = znimg.darr


    # Determine number of tiles in each dimension
    img_shape = img_data.shape
    step = tuple(s - overlap for s in tile_shape)

    grid_shape = tuple(
        max(((img_shape[dim] - tile_shape[dim]) // step[dim]) + 1, 1)
        for dim in range(3)
    )

    print(f'img_shape {img_shape}, step: {step}, grid_shape: {grid_shape}')
    # Create tiles
    tiles = []
    translations = []
    for z in range(grid_shape[0]):
        for y in range(grid_shape[1]):
            for x in range(grid_shape[2]):
                # Extract tile
                z_start, y_start, x_start = z * step[0], y * step[1], x * step[2]
                tile = img_data[z_start:z_start+tile_shape[0],
                                y_start:y_start+tile_shape[1],
                                x_start:x_start+tile_shape[2]]

                # Add to list
                tiles.append(tile)

                # Add random offset
                offset = np.random.uniform(-5, 5, size=3)  # Random 3D offsets
                translations.append((z_start, y_start, x_start) + offset)

    print(tiles)
    # Convert to a Dask array
    print(tile_shape)
    tiles = da.concatenate([tile.rechunk(chunks=final_chunks) for tile in tiles])
    translations = np.array(translations)

    znimg.darr = tiles

    return znimg, translations




if __name__ == '__main__':
    test_znimg, test_translations = create_test_dataset()
    print(test_znimg)

    print(test_translations.shape)
    test_znimg.to_ome_zarr(out_ome_zarr)
    np.save(out_translations_npy,test_translations)
    
