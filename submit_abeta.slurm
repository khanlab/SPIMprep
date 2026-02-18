#!/bin/bash
#SBATCH --job-name=spimprep_abeta
#SBATCH --output=logs/abeta_%A_%a.out
#SBATCH --error=logs/abeta_%A_%a.err
#SBATCH --array=1-4%2
#SBATCH --time=8:00:00
#SBATCH --mem=250G
#SBATCH --cpus-per-task=32
#SBATCH --export=ALL

SPIMPREP_DIR="/nfs/khan/datasets/selma_all/SPIMprep"
RAW_DIR="/nfs/khan/datasets/selma_all/Ab_plaques/raw"
BIDS_DIR="/nfs/khan/datasets/selma_all/Ab_plaques/bids"
WORK_DIR="/tmp/spimprep_abeta_${SLURM_ARRAY_TASK_ID}"

SAMPLE_NUM=${SLURM_ARRAY_TASK_ID}
SUBJECT=$(printf "%02d" ${SAMPLE_NUM})

mkdir -p "${WORK_DIR}"

cd "${SPIMPREP_DIR}"

pixi run ./run.py \
  --input-path "${RAW_DIR}/sample${SAMPLE_NUM}" \
  --subject "${SUBJECT}" \
  --sample brain \
  --acq prestitched \
  --stains abeta \
  --resolution-preset selma_abeta \
  --output-bids-dir "${BIDS_DIR}" \
  --work-dir "${WORK_DIR}" \
  --cores ${SLURM_CPUS_PER_TASK}

rm -rf "${WORK_DIR}"