# Snakemake workflow: `SPIMprep`

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/khanlab/SPIMprep/workflows/Tests/badge.svg?branch=main)](https://github.com/khanlab/SPIMprep/actions?query=branch%3Amain+workflow%3ATests)

A Snakemake workflow for pre-processing single plane illumination microscopy (SPIM, aka lightsheet microscopy)

Takes datasets produced by the Biotec UltraMicroscope Blaze (as a folder of tif images), and outputs a BIDS-like dataset, with a multi-channel multi-scale ome zarr file along with downsampled nifti images. 

## Requirements

 - Linux system with Singularity/Apptainer installed 
   - (Note: container will be automatically pulled when you run the workflow)
 - Python >= 3.11
 - Lightsheet data as OME tif files in a folder, or in a tar file (see tif filename pattern in `config.yml`)

## Usage


1. Clone this repository to the folder you want to run the workflow in
```
git clone https://github.com/khanlab/spimprep
```

2. Create and activate a virtual environment, then install dependencies with:

pip install .
```
Note: to make a venv on the CBS server use:
```
python3.11 -m venv venv
source venv/bin/activate
```

3. Update the `config/datasets.tsv` spreadsheet to point to your dataset(s). Each dataset's tif files should be in it's own folder or tar file, with no other tif files. Enter the path to each dataset in the `dataset_path` column. The first three columns identify the subject, sample, acquisition, which become part of the resulting filenames (BIDS naming). The `stain_0` and `stain_1` identify what stains were used for each channel. Use `autof` to indicate the autofluorescence channel, as this is used for registration.

4. The `config/config.yml` can be edited to customize any workflow parameters. The most important ones are the `root` and `work` variables. The `root` path is where the results will end up, by default this is a subfolder called results. The `work` path is where any intermediate scratch files are produced. By default the files in `work` are deleted after they are no longer needed in the workflow, unless you use the `--notemp` command-line option. The workflow writes a large number of small files in parallel to the `work` folder, so for optimum performance this should be a fast local disk, and not a networked file system (i.e. shared disk).  

Note: you can use environment variables when specifying `root` or `work`, e.g. so `work: '$SLURM_TMPDIR` can be used on HPC servers. 
5. Perform a dry-run to make sure the workflow is configured properly. This will only print what the workflow will run, and will not run anything.
```
snakemake -np 
```

6.  To run the workflow, parallelizing on all cores, using Singularity (aka Apptainer) for dependencies, use:
```
snakemake -c all --sdm apptainer 
```

Note: if you run the workflow on a system with large memory, you will need to set the heap size for the stitching and fusion rules. This can be done with e.g.: `--set-resources bigstitcher:mem_mb=60000 fuse_dataset:mem_mb=100000`

7. If you want to run the workflow using a batch job submission server, please see the executor plugins here: https://snakemake.github.io/snakemake-plugin-catalog/



### WIP:  The usage of this workflow is described in the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/?usage=khanlab%2Fspimprep).


