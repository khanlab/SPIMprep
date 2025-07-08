#!/bin/bash
#SBATCH --nodes=1

#eval "$(pixi shell-hook)"

# Optionally print some env info
echo "Environment activated with pixi."
echo "Python version: $(python --version)"
echo "Using RUN_PY_PATH: $RUN_PY_PATH"

"$RUN_PY_PATH" \
  --output-bids-dir $OUTPUT_BIDS_DIR \
  --work-dir $WORK_DIR \
  --stains $STAINS \
  --subject $SUBJECT \
  --acq $ACQ \
  --input-path $INPUT_PATH \
  -c all \
  --use-conda \
  --conda-prefix $CONDA_PREFIX


