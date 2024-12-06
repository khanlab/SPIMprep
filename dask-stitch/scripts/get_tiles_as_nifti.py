import nibabel as nib
from zarrnii import ZarrNii
from pathlib import Path


out_dir = Path(snakemake.output.tiles_dir)
out_dir.mkdir(exist_ok=True, parents=True)

for tile in range(snakemake.params.n_tiles):
    print(f'reading tile {tile} and writing to nifti')
    ZarrNii.from_path(snakemake.input.ome_zarr,channels=[tile]).to_nifti(out_dir / f'tile_{tile:02d}.nii')



