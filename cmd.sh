#!/bin/bash
unset SNAKEMAKE_PROFILE
snakemake --jobs 100 --executor googlebatch --googlebatch-region us-east1 --googlebatch-project t-system-193821 --default-storage-provider gcs --default-storage-prefix gs://khanlab-bucket/spimprep/test_run_8000milli --storage-gcs-project t-system-193821  --sdm apptainer conda -p  --googlebatch-image-project t-system-193821  --googlebatch-boot-disk-image projects/t-system-193821/global/images/custom-batch-image-apptainer --googlebatch-boot-disk-gb 150 --googlebatch-machine-type c2-standard-8 --googlebatch-cpu-milli 8000 --googlebatch-memory 16000 $@
