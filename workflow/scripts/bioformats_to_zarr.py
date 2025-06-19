from pathlib import Path
import subprocess
import shlex

# Folder containing the .ome.tif files
#input_dir = Path(snakemake.input.ome_dir)
input_dir = Path(snakemake.params.ome_dir)
output_base = Path(snakemake.output.tiles_dir)
output_base.mkdir(parents=True, exist_ok=True)

# Loop through all files
for fname in sorted(input_dir.iterdir()):
    if fname.name.endswith("C00.ome.tif"):
        input_path = fname
        
        # Strip extension for output subdir
        base_name = fname.name.rsplit(".ome.tif", 1)[0]
        output_bf2raw = output_base / f'{base_name}.bf2raw'

        # Construct command
        cmd = f'bioformats2raw "{input_path}" "{output_bf2raw}" -h {snakemake.params.tile_height} -w {snakemake.params.tile_width} -p -r 1 --max-workers {snakemake.threads}'
        print(f"Running: {cmd}")
        subprocess.run(shlex.split(cmd), check=True)

