import nibabel as nib
import numpy as np
from scipy.ndimage import affine_transform
from templateflow import api as tflow
import nibabel as nib
import dask.array as da
import math
from zarrnii import ZarrNii

def create_test_dataset_single(tile_index, template="MNI152NLin2009cAsym", res=2, grid_shape=(3, 4), overlap=32, random_seed=42, final_chunks=(1,32, 32, 1)):
    """
    Create a low-resolution test dataset for tile-based stitching.

    Parameters:
    - tile_index: the index of the tile to create
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
    img_data = nib.load(template_path).get_fdata()

    # Original image shape
    img_shape = np.array(img_data.shape)  # (Z, Y, X)

    # Keep Z dimension intact, calculate X and Y tile sizes
    x_dim, y_dim, z_dim = img_shape
    x_tile_size = math.ceil((x_dim + overlap * (grid_shape[0] - 1)) / grid_shape[0])
    y_tile_size = math.ceil((y_dim + overlap * (grid_shape[1] - 1)) / grid_shape[1])

    # Calculate required padding to ensure X and Y dimensions are divisible by grid shape
    padded_x = x_tile_size * grid_shape[0] - overlap * (grid_shape[0] - 1)
    padded_y = y_tile_size * grid_shape[1] - overlap * (grid_shape[1] - 1)

    padding = (
        (0, int(max(padded_x - x_dim, 0))),
        (0, int(max(padded_y - y_dim, 0))),
        (0, 0),  # No padding in Z
    )


    # Pad image if needed
    if any(p[1] > 0 for p in padding):
        img_data = np.pad(img_data, padding, mode="constant", constant_values=0)

    # Create tiles

    x,y = np.unravel_index(tile_index,grid_shape)


    # Calculate tile start indices
    x_start = x * (x_tile_size - overlap)
    y_start = y * (y_tile_size - overlap)


    # TODO: Simulate error by applying a transformation to the image before 

    # initially lets just do a random jitter:
    offset = np.random.uniform(0, 0, size=(grid_shape[0],grid_shape[1],3))  # Random 3D offsets for each tile

    xfm_img_data = affine_transform(img_data,matrix=np.eye(3,3),offset=offset[x,y,:],order=1)

    # Extract tile
    tile = xfm_img_data[x_start:x_start + x_tile_size, y_start:y_start + y_tile_size, :]



    translation = ((x_start, y_start, 0))

    

    tile_shape = (1,x_tile_size, y_tile_size, z_dim)


    #save tiling coordinate translation into vox2ras (not the random jitter)
    vox2ras = np.eye(4)
    vox2ras[:3,3] = np.array(translation)

    #this next bit is hacky and in future should be incorporated into ZarrNii.from_darr() 
    reorder_xfm = np.eye(4)
    reorder_xfm[:3, :3] = np.flip(
        reorder_xfm[:3, :3], axis=0
    )  # reorders z-y-x to x-y-z and vice versa

    flip_xfm = np.diag((-1,-1,-1,1))
    
    vox2ras = flip_xfm @ reorder_xfm @ vox2ras
        

    # Save back into ZarrNii object
    darr = da.from_array(tile.reshape(tile_shape),chunks=final_chunks)
   
    
    znimg = ZarrNii.from_darr(darr,vox2ras=vox2ras)

    return znimg




test_znimg = create_test_dataset_single(tile_index=snakemake.params.tile_index,
                                                    grid_shape=snakemake.params.grid_shape)
test_znimg.to_ome_zarr(snakemake.output.ome_zarr)
test_znimg.to_nifti(snakemake.output.nifti)



print("LOOP0")
print(test_znimg)

test_znimg2 = ZarrNii.from_path(snakemake.output.ome_zarr)
print("LOOP1")
print(test_znimg2) #we want this to be identical to test_znimg

#the zarr we read in has a different vox2ras - the translation is negated and flipped
updated_affine = test_znimg2.vox2ras.affine.copy()
print(f'original affine: {updated_affine}')

#HACK FIX:
#affine has translation negated and flipped, for whatever reason..  need to fix this in zarrnii
updated_affine[:3,3] = -1 * np.flip(updated_affine[:3,3])

test_znimg2.vox2ras.affine = updated_affine
test_znimg2.ras2vox.affine = np.linalg.inv(updated_affine)


test_znimg2.to_nifti(snakemake.output.nifti_fromzarr)

"""

test_znimg2.to_ome_zarr('mytest_writeafterread.ome.zarr')
test_znimg3 = ZarrNii.from_path('mytest_writeafterread.ome.zarr')
print("LOOP2")
print(test_znimg3)


test_znimg3.to_ome_zarr('mytest_writeafterread2.ome.zarr')
test_znimg4 = ZarrNii.from_path('mytest_writeafterread2.ome.zarr')
print("LOOP3")
print(test_znimg4)
"""
