import nibabel as nib
from zarrnii import ZarrNii
from pathlib import Path


znimg_example_tile= ZarrNii.from_path(snakemake.input.ome_zarr)
print(znimg_full.darr.shape)

out_dir = Path(snakemake.output.tiles_dir)
out_dir.mkdir(exist_ok=True, parents=True)

for tile in range(znimg_full.darr.shape[0]):
    print(f'reading tile {tile} and writing to nifti')
    ZarrNii.from_path(snakemake.input.ome_zarr,channels=[tile]).to_nifti(out_dir / f'tile_{tile:02d}.nii')



