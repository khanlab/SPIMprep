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

pixi run python << 'PYEOF'
import os
import sys
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

# Compute downsample factors to get close to 4 µm (0.004 mm) isotropic
factors = []
for dim in ["z", "y", "x"]:
    ratio = target_mm / scale[dim]
    if ratio >= 1.5:
        # Use nearest power of 2
        factor = max(1, 2 ** int(np.round(np.log2(ratio))))
    else:
        factor = 1
    factors.append(factor)

print(f"Downsample factors (z, y, x): {factors}")

if any(f > 1 for f in factors):
    znimg = znimg.downsample(factors=factors, spatial_dims=["z", "y", "x"])

final_scale = znimg.scale
print(f"Final scales (mm): z={final_scale['z']:.6f}, y={final_scale['y']:.6f}, x={final_scale['x']:.6f}")
print(f"Final scales (µm): z={final_scale['z']*1000:.2f}, y={final_scale['y']*1000:.2f}, x={final_scale['x']*1000:.2f}")
print(f"Final shape: {znimg.darr.shape}")

print(f"Writing to: {zarr_out}")
with ProgressBar():
    znimg.to_ome_zarr(zarr_out, max_layer=4, backend="ome-zarr-py", scaling_method="local_mean")

print("Done!")
PYEOF

echo "Finished ${SUB_ID}"
