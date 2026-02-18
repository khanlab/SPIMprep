# %% Test if the prestitched codepath imports will work
import zarr
print(f"zarr version: {zarr.__version__}")

# These are what tif_stacks_to_ome_zarr.py needs:
try:
    from ome_zarr.writer import write_image
    from ome_zarr.scale import Scaler
    from ome_zarr.io import parse_url
    print("ome_zarr imports: OK")
except ImportError as e:
    print(f"ome_zarr imports FAILED: {e}")

# Test the v2 store APIs the script uses:
try:
    store = zarr.DirectoryStore("/tmp/test_zarr", dimension_separator="/")
    print("zarr.DirectoryStore: OK")
except AttributeError:
    print("zarr.DirectoryStore: MISSING (zarr v3 breaking change)")

# Also test zarrnii (used by ome_zarr_to_nii for nifti export):
try:
    from zarrnii import ZarrNii
    print("zarrnii imports: OK")
except ImportError as e:
    print(f"zarrnii imports FAILED: {e}")
