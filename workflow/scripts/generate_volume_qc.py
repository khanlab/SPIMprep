import json
import shutil
from pathlib import Path
from distutils.dir_util import copy_tree
from zarrnii import ZarrNii
import math

# directory containing the volume rendering files
resource_dir = Path(snakemake.output.resources)
# where html file should be written
html_dest = snakemake.output.html

# inputted ome-zarr path
ome_data = snakemake.input.ome

# move volume renderer into the subjects directory
copy_tree("qc/resources/volViewer", str(resource_dir))
shutil.move(resource_dir / "volRender.html", html_dest)

# Get most downsampled ome-zarr image
ds_z = ZarrNii.from_path(ome_data,level=5, channels=[0,1])
z_length = ds_z.darr.shape[1]

# downsample it so it has at most 100 slices and ast least 50 slices in z-direction
if(z_length>50):
    downsample_factor = math.floor(z_length/50)
else:
    downsample_factor = 1
ds_z.downsample(along_z=downsample_factor)

# Write it to a JSON for js script to read
with open(resource_dir / "volumeData.json", 'w') as f:
    json_data = json.dumps(ds_z.darr.compute().tolist())
    f.write(json_data)




