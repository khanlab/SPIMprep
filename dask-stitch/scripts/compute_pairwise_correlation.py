import nibabel as nib
import numpy as np
from zarrnii import ZarrNii, AffineTransform
from scipy.fft import fftn, ifftn
from scipy.ndimage import map_coordinates
import os
from skimage.registration import phase_cross_correlation

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
    print(f'Performing phase_correlation on pair {pair_index}')

    # Compute shifts using skimage's phase_cross_correlation
    shifts, error, diffphase = phase_cross_correlation(img1, img2, upsample_factor=10, space='real', disambiguate=True)

    print(f"Shifts: {shifts}, Error: {error}, Phase Difference: {diffphase}")

    # Optionally save a diagnostic correlation map
    if diagnostics_dir and pair_index is not None:
        # Generate a mock correlation map as skimage does not provide it directly
        correlation_map = np.roll(np.roll(np.roll(img1, int(shifts[0]), axis=0), int(shifts[1]), axis=1), int(shifts[2]), axis=2)
        np.save(os.path.join(diagnostics_dir, f"correlation_map_pair_{pair_index}.npy"), correlation_map)
        nib.Nifti1Image(correlation_map, affine=np.eye(4)).to_filename(
            f'{diagnostics_dir}/correlation_map_pair_{pair_index}.nii'
        )

    return shifts


def crop_to_bounding_box_pair(img1, affine1, img2, affine2, bbox_min, bbox_max):
    """
    Crop two images to the specified bounding box in physical space, ensuring identical sizes.

    Parameters:
    - img1 (np.ndarray): First input image.
    - affine1 (np.ndarray): 4x4 affine matrix mapping voxel to physical space for the first image.
    - img2 (np.ndarray): Second input image.
    - affine2 (np.ndarray): 4x4 affine matrix mapping voxel to physical space for the second image.
    - bbox_min (np.ndarray): Minimum physical coordinates of the bounding box.
    - bbox_max (np.ndarray): Maximum physical coordinates of the bounding box.
    - output_shape (tuple, optional): Desired shape of the output images.

    Returns:
    - tuple of np.ndarray: Cropped images (cropped_img1, cropped_img2).
    """
    def get_voxel_bounds(img, affine, bbox_min, bbox_max):
        """Transform physical bounding box to voxel coordinates and clip to image bounds."""
        affine_inverse = AffineTransform.from_array(affine,invert=True)

        # Transform physical to voxel space
        voxel_min = np.floor(affine_inverse @ bbox_min).astype(int)
        voxel_max = np.ceil(affine_inverse @ bbox_max).astype(int)

        # Clip to image bounds
        voxel_min = np.clip(voxel_min, 0, np.array(img.shape) - 1)
        voxel_max = np.clip(voxel_max, 0, np.array(img.shape) - 1)

        return voxel_min, voxel_max

    # Get voxel bounds for both images
    voxel_min1, voxel_max1 = get_voxel_bounds(img1, affine1, bbox_min, bbox_max)
    voxel_min2, voxel_max2 = get_voxel_bounds(img2, affine2, bbox_min, bbox_max)


    # Ensure the cropped regions have the same physical size
    size1 = voxel_max1 - voxel_min1
    size2 = voxel_max2 - voxel_min2

    # Adjust voxel_max2 to match the physical size of the region in img1
    size = np.minimum(size1, size2)
    voxel_max1 = voxel_min1 + size
    voxel_max2 = voxel_min2 + size

    # Crop the images
    cropped_img1 = img1[
        voxel_min1[0]:voxel_max1[0] + 1,
        voxel_min1[1]:voxel_max1[1] + 1,
        voxel_min1[2]:voxel_max1[2] + 1,
    ]
    cropped_img2 = img2[
        voxel_min2[0]:voxel_max2[0] + 1,
        voxel_min2[1]:voxel_max2[1] + 1,
        voxel_min2[2]:voxel_max2[2] + 1,
    ]

    return cropped_img1, cropped_img2



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
        np.array([img_shape[0], 0, 0, 1]),
        np.array([0, img_shape[1], 0, 1]),
        np.array([img_shape[0], img_shape[1], 0, 1]),
        np.array([0, 0, img_shape[2], 1]),
        np.array([img_shape[0], 0, img_shape[2], 1]),
        np.array([0, img_shape[1], img_shape[2], 1]),
        np.array([img_shape[0], img_shape[1], img_shape[2], 1]),
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
        znimg1 = ZarrNii.from_ome_zarr(ome_zarr_paths[i])
        znimg2 = ZarrNii.from_ome_zarr(ome_zarr_paths[j])

        img1 = znimg1.darr.squeeze().compute()
        img2 = znimg2.darr.squeeze().compute()

        #TODO: should make this a class member function, or just return the reordered affine
        affine1 = np.eye(4,4)
        affine1[:3, :3] = np.diag(znimg1.get_zooms(axes_order=znimg1.axes_order))  # Set the zooms (scaling factors) along the diagonal
        affine1[:3, 3] = znimg1.get_origin(axes_order=znimg1.axes_order)  # Set the translation (origin)

        affine2 = np.eye(4,4)
        affine2[:3, :3] = np.diag(znimg2.get_zooms(axes_order=znimg2.axes_order))  # Set the zooms (scaling factors) along the diagonal
        affine2[:3, 3] = znimg2.get_origin(axes_order=znimg2.axes_order)  # Set the translation (origin)


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
        resampled_img1, resampled_img2 = crop_to_bounding_box_pair(img1, affine1, img2, affine2, bbox_min, bbox_max)
      #  resampled_img1 = resample_to_bounding_box(img1, bbox_min, bbox_max, output_shape)
      #  resampled_img2 = resample_to_bounding_box(img2, bbox_min, bbox_max, output_shape)

        # Save resampled images for diagnostics
        if diagnostics_dir:
            nib.Nifti1Image(resampled_img1,affine=np.eye(4,4)).to_filename(f'{diagnostics_dir}/resampled_img1_pair_{pair_index}.nii')
            nib.Nifti1Image(resampled_img2,affine=np.eye(4,4)).to_filename(f'{diagnostics_dir}/resampled_img2_pair_{pair_index}.nii')

        # Compute phase correlation on the resampled overlapping region
        offset = phase_correlation(resampled_img1, resampled_img2, diagnostics_dir=diagnostics_dir, pair_index=pair_index)
        if znimg1.axes_order == 'ZYX':
            offset = offset[::-1]
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


