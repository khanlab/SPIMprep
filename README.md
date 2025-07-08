# Snakemake workflow: `SPIMprep`

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/khanlab/SPIMprep/workflows/Tests/badge.svg?branch=main)](https://github.com/khanlab/SPIMprep/actions?query=branch%3Amain+workflow%3ATests)

<!--intro-start-->

A Snakemake workflow for pre-processing single plane illumination microscopy (SPIM, aka lightsheet microscopy).

Takes TIF images (tiled or prestitched) and outputs a validated BIDS Microscopy dataset, with a multi-channel multi-scale OME-Zarr file for each scan, along with downsampled nifti images (in a derivatives folder). 

## Supported inputs:

SPIMprep supports a range of inputs, with the type of acquisition specified by including
the short-hand name (in bold below) as a substring in the acquisition tag.
  - **`blaze`**:  Raw Ultramicroscope Blaze OME TIFF files, either as 2D or 3D TIFF files
  - **`prestitched`**: Prestitched images, as a stack of 2D TIF files (e.g. from LifeCanvas)
  - **`imaris`**: Prestitched into a single Imaris (.ims) file.


## Requirements

 - Linux system with Singularity/Apptainer installed 
    - (Note: container will be automatically pulled when you run the workflow)
 - Python >= 3.11
 - Lightsheet data:

## Usage


1. Clone this repository to the folder you want to run the workflow in
```
git clone https://github.com/khanlab/spimprep
```

2. Create and activate a virtual environment, then install dependencies with:
```
pip install .
```
Note: to make a venv on the CBS server use:
```
python3.11 -m venv venv
source venv/bin/activate
```

3. Update the `config/samples.tsv` spreadsheet to point to your sample(s). Each sample's tif files should be in it's own folder or tar file, with no other tif files. Enter the path to each sample in the `sample_path` column. The first three columns identify the subject, sample, acquisition, which become part of the resulting filenames (BIDS naming). The `stain_0` and `stain_1` identify what stains were used for each channel. Use `autof` to indicate the autofluorescence channel. If you have a different number of stains you can add or remove these columns. If your samples have different numbers of stains, you can leave values blank or use `n/a` to indicate that a sample does not have a particular stain. 

Note: The acquisition value must contain either `blaze` or `prestitched`, and defines which workflow will be used. E.g. for LifeCanvas data that is already stitched, you need to include `prestitched` in the acquisition flag. 

**New:** Writing output directly to cloud storage is now supported; enable this by using `s3://` or `gcs://` in the `root` variable, to point to a bucket you have write access to. 

5. The `config/config.yml` can be edited to customize any workflow parameters. The most important ones are the `root` and `work` variables. The `root` path is where the results will end up, by default this is a subfolder called `bids`. The `work` path is where any intermediate scratch files are produced. By default the files in `work` are deleted after they are no longer needed in the workflow, unless you use the `--notemp` command-line option. The workflow writes a large number of small files in parallel to the `work` folder, so for optimum performance this should be a fast local disk, and not a networked file system (i.e. shared disk).  

Note: you can use environment variables when specifying `root` or `work`, e.g. so `work: '$SLURM_TMPDIR` can be used on HPC servers. 

5. Go to the SPIMprep folder and perform a dry-run to make sure the workflow is configured properly. This will only print what the workflow will run, and will not run anything.
```
snakemake -np 
```

6.  To run the workflow, parallelizing on all cores, using Singularity (aka Apptainer) for dependencies, use:
```
snakemake -c all --sdm apptainer 
```
or for snakemake<8.0, use:
```
snakemake -c all --use-singularity 
```

Note: if you run the workflow on a system with large memory, you will need to set the heap size for the stitching and fusion rules. This can be done with e.g.: `--set-resources bigstitcher_stitching:mem_mb=60000 bigstitcher_fusion:mem_mb=100000`

7. If you want to run the workflow using a batch job submission server, please see the executor plugins here: https://snakemake.github.io/snakemake-plugin-catalog/

<!--intro-end-->



## To incorporate later (new CLI):

```
git clone https://github.com/khanlab/spimprep
git checkout multiview-stitcher
pixi install
pixi shell-hook > activate
chmod a+x ./activate
source ./activate
./run.py  --input-path /cifs/khan_new/datasets/lightsheet_example/ --subject M4A1Te3 --stains Abeta PI Lectin --output-bids-dir /tmp/out_spimprep --work-dir /tmp   --cores all
