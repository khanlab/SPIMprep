
import numpy as np
from zarrnii import ZarrNii, AffineTransform

def assign_translations(ome_zarr_paths, translations, output_niftis, output_ome_zarrs):
    """
    Update the OME-Zarr datasets with optimized translations by modifying the vox2ras.

    Parameters:
    - ome_zarr_paths (list of str): List of input OME-Zarr dataset paths.
    - translations (np.ndarray): Optimized translations (T, 3).
    - output_niftis (list of str): List of paths to save the updated nifti datasets.
    - output_ome_zarrs (list of str): List of paths to save the updated OME-Zarr datasets.
    """
    if len(ome_zarr_paths) != len(translations):
        raise ValueError("Number of OME-Zarr paths must match the number of optimized translations.")

    for i, (path, translation, out_nii, out_zarr) in enumerate(zip(ome_zarr_paths, translations, output_niftis, output_ome_zarrs)):
        print(path)
        print(translation)
        # Load the OME-Zarr dataset
        znimg = ZarrNii.from_ome_zarr(path)

        affine = znimg.affine
   
        if znimg.axes_order == 'ZYX':
            affine[:3, 3] += translation[::-1]  # Add the optimized translation
        else:
            affine[:3, 3] += translation  # Add the optimized translation

        znimg.affine = affine

        # Save the updated ZarrNii
        znimg.to_nifti(out_nii)
        znimg.to_ome_zarr(out_zarr)

        print(f"Updated translations saved to: {out_nii} and {out_zarr}")


# Example usage
ome_zarr_paths = snakemake.input.ome_zarr  # List of input OME-Zarr paths
translations = np.loadtxt(snakemake.input.translations, dtype=float)  # Optimized translations
output_niftis = snakemake.output.niftis  # List of output nifti paths
output_ome_zarrs = snakemake.output.ome_zarrs  # List of output nifti paths

# Assign translations
assign_translations(ome_zarr_paths, translations, output_niftis, output_ome_zarrs)

