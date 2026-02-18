#!/bin/bash
#SBATCH --job-name=spimprep_vessel
#SBATCH --output=logs/vessel_%A_%a.out
#SBATCH --error=logs/vessel_%A_%a.err
#SBATCH --array=1-18%2
#SBATCH --time=8:00:00
#SBATCH --mem=250G
#SBATCH --cpus-per-task=32
#SBATCH --export=ALL

SPIMPREP_DIR="/nfs/khan/datasets/selma_all/SPIMprep"
RAW_DIR="/nfs/khan/datasets/selma_all/vessel/raw"
BIDS_DIR="/nfs/khan/datasets/selma_all/vessel/bids"
WORK_DIR="/tmp/spimprep_vessel_${SLURM_ARRAY_TASK_ID}"

SAMPLE_NUM=${SLURM_ARRAY_TASK_ID}
SUBJECT=$(printf "%02d" ${SAMPLE_NUM})

mkdir -p "${WORK_DIR}"

cd "${SPIMPREP_DIR}"

pixi run ./run.py \
  --input-path "${RAW_DIR}/sample${SAMPLE_NUM}" \
  --subject "${SUBJECT}" \
  --sample brain \
  --acq prestitched \
  --stains vessel \
  --resolution-preset selma_blood_vessel \
  --output-bids-dir "${BIDS_DIR}" \
  --work-dir "${WORK_DIR}" \
  --cores ${SLURM_CPUS_PER_TASK}

rm -rf "${WORK_DIR}"
