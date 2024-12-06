
import numpy as np
from zarrnii import ZarrNii

def assign_translations(ome_zarr_paths, optimized_translations, output_paths):
    """
    Update the OME-Zarr datasets with optimized translations by modifying the vox2ras.

    Parameters:
    - ome_zarr_paths (list of str): List of input OME-Zarr dataset paths.
    - optimized_translations (np.ndarray): Optimized translations (T, 3).
    - output_paths (list of str): List of paths to save the updated OME-Zarr datasets.
    """
    if len(ome_zarr_paths) != len(optimized_translations):
        raise ValueError("Number of OME-Zarr paths must match the number of optimized translations.")

    for i, (path, translation, output_path) in enumerate(zip(ome_zarr_paths, optimized_translations, output_paths)):
        # Load the OME-Zarr dataset
        znimg = ZarrNii.from_path(path)

        # Update the affine matrix
        updated_affine = znimg.vox2ras.affine.copy()
        updated_affine[:3, 3] += translation  # Add the optimized translation

        # Create a new ZarrNii object with the updated affine
        updated_znimg = ZarrNii.from_darr(
            darr=znimg.darr,
            vox2ras=updated_affine,
            axes_nifti=znimg.axes_nifti  # Preserve existing NIFTI axis convention
        )

        # Save the updated ZarrNii
        updated_znimg.to_nifti(output_path)
        print(f"Updated translations saved to: {output_path}")


# Example usage
ome_zarr_paths = snakemake.input.ome_zarr  # List of input OME-Zarr paths
optimized_translations = np.loadtxt(snakemake.input.optimized_translations, dtype=float)  # Optimized translations
output_paths = snakemake.output.nifti  # List of output nifti paths

# Assign translations
assign_translations(ome_zarr_paths, optimized_translations, output_paths)

