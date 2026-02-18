#!/bin/bash
#SBATCH --job-name=ds4um_abeta
#SBATCH --output=logs/ds4um_abeta_%A_%a.out
#SBATCH --error=logs/ds4um_abeta_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --array=0-3

# ---- Setup ----
BASE_DIR="/nfs/khan/datasets/selma_all"
DATASET="Ab_plaques"
BIDS_DIR="${BASE_DIR}/${DATASET}/bids"
OUT_DIR="${BASE_DIR}/${DATASET}/bids/derivatives/4um_isotropic"

mkdir -p logs

# Build array of subject directories
SUBJECTS=($(ls -d ${BIDS_DIR}/sub-* | sort))
SUB=${SUBJECTS[$SLURM_ARRAY_TASK_ID]}
SUB_ID=$(basename ${SUB})

echo "Processing ${SUB_ID} (task ${SLURM_ARRAY_TASK_ID})"

# Find the OME-Zarr for this subject
ZARR_IN=$(ls -d ${SUB}/micr/*.ome.zarr 2>/dev/null | head -1)

if [ -z "${ZARR_IN}" ]; then
    echo "No OME-Zarr found for ${SUB_ID}, skipping."
    exit 0
fi

ZARR_BASENAME=$(basename ${ZARR_IN})
OUT_SUB_DIR="${OUT_DIR}/${SUB_ID}/micr"
ZARR_OUT="${OUT_SUB_DIR}/${ZARR_BASENAME}"

if [ -d "${ZARR_OUT}" ]; then
    echo "Output already exists: ${ZARR_OUT}, skipping."
    exit 0
fi

mkdir -p ${OUT_SUB_DIR}

echo "Input:  ${ZARR_IN}"
echo "Output: ${ZARR_OUT}"

export ZARR_IN ZARR_OUT

python3 << 'PYEOF'
import os
import numpy as np
from zarrnii import ZarrNii
from dask.diagnostics import ProgressBar

zarr_in = os.environ["ZARR_IN"]
zarr_out = os.environ["ZARR_OUT"]
target_mm = 0.004  # 4 µm in mm

print(f"Loading level 0: {zarr_in}")
znimg = ZarrNii.from_ome_zarr(zarr_in, level=0)

scale = znimg.scale  # in mm
print(f"Original scales (mm): z={scale['z']:.6f}, y={scale['y']:.6f}, x={scale['x']:.6f}")
print(f"Original scales (µm): z={scale['z']*1000:.2f}, y={scale['y']*1000:.2f}, x={scale['x']*1000:.2f}")
print(f"Original shape: {znimg.darr.shape}")
print(f"Dims: {znimg.dims}")

# Compute exact zoom factors to reach 4 µm isotropic
# zoom < 1 = downsampling, zoom > 1 = upsampling
zoom_factors = {dim: scale[dim] / target_mm for dim in ["z", "y", "x"]}
print(f"Zoom factors: z={zoom_factors['z']:.4f}, y={zoom_factors['y']:.4f}, x={zoom_factors['x']:.4f}")

# Build target shape: keep non-spatial dims, rescale spatial dims
target_shape = []
for dim in znimg.dims:
    dl = dim.lower()
    dim_idx = list(d.lower() for d in znimg.dims).index(dl)
    orig_size = znimg.darr.shape[dim_idx]
    if dl in zoom_factors:
        target_shape.append(int(np.round(orig_size * zoom_factors[dl])))
    else:
        target_shape.append(orig_size)

target_shape = tuple(target_shape)
print(f"Target shape: {target_shape}")

# upsample with to_shape uses scipy.ndimage.zoom internally — works for
# both upsampling and downsampling despite the method name
znimg_resampled = znimg.upsample(to_shape=target_shape)

final_scale = znimg_resampled.scale
print(f"Final scales (mm): z={final_scale['z']:.6f}, y={final_scale['y']:.6f}, x={final_scale['x']:.6f}")
print(f"Final scales (µm): z={final_scale['z']*1000:.2f}, y={final_scale['y']*1000:.2f}, x={final_scale['x']*1000:.2f}")
print(f"Final shape: {znimg_resampled.darr.shape}")

print(f"Writing to: {zarr_out}")
with ProgressBar():
    znimg_resampled.to_ome_zarr(zarr_out, max_layer=4, backend="ome-zarr-py", scaling_method="local_mean")

print("Done!")
PYEOF

echo "Finished ${SUB_ID}"
