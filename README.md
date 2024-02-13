# Snakemake workflow: `SPIMprep`

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/<owner>/<repo>/workflows/Tests/badge.svg?branch=main)](https://github.com/<owner>/<repo>/actions?query=branch%3Amain+workflow%3ATests)

A Snakemake workflow for pre-processing single plane illumination microscopy (SPIM, aka lightsheet microscopy)

Takes datasets produced by the Biotec UltraMicroscope Blaze (as a folder of tif images), and outputs a BIDS-like dataset, with multiscale ome zarr files along with downsampled nifti images. 

Also performs simple affine registration to Allen Brain Atlas v3 (mouse atlas), and resamples labels to multiscale ome zarr as well. (WIP)

## Requirements

 - Linux system with Singularity/Apptainer installed.
 - Python >= 3.11
 - Lightsheet data in a folder, or in a tar file (see tif filename pattern in `config.yml`)

## Usage


Clone this repository to the folder you want to run the workflow in, update the `config/datasets.tsv` to point to your dataset(s), and optionally update the `config/config.yml` with customized options. 

Install dependencies with:
```
pip install .
```

Then, to do a dry-run use:
```
snakemake -np 
```

To run on all cores use:
```
snakemake -c all --sdm apptainer
```

### WIP:  The usage of this workflow is described in the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/?usage=khanlab%2Fspimprep).



## TODO:

 - [ ] any changes to post to Snakemake Workflow Catalog, and verify usage there
