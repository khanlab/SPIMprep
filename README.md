# `SPIMprep`

<!--intro-start-->

A Snakemake+SnakeBIDS workflow for pre-processing single plane illumination microscopy (SPIM, aka lightsheet microscopy).

Takes TIF images (tiled or prestitched) and outputs a validated BIDS Microscopy dataset, with a multi-channel multi-scale OME-Zarr file for each scan, along with downsampled nifti images (in a derivatives folder). 

## Supported inputs:

SPIMprep supports a range of inputs, with the type of acquisition specified by including
the short-hand name (in bold below) as a substring in the acquisition tag.
  - **`blaze`**:  Raw Ultramicroscope Blaze OME TIFF files, either as 2D or 3D TIFF files
  - **`prestitched`**: Prestitched images, as a stack of 2D TIF files (e.g. from LifeCanvas)
  - **`imaris`**: Prestitched into a single Imaris (.ims) file.


## Requirements

 - Pixi package manager (https://pixi.sh/latest/)

## Installation



0. If pixi is not installed, install it as described here: https://pixi.sh/latest/installation/

1. Clone this repository and cd to it:
```
git clone https://github.com/khanlab/spimprep SPIMprep
cd SPIMprep
```

2. Install any dependencies using pixi:
```
pixi install
```

3. Activate the created environment using pixi, or create an activate script (if `pixi shell` does not work):
```
pixi shell
```
or
```
pixi shell-hook > activate
source ./activate
```

## Usage

Run the workflow, using the `run.py` script:
```
usage: run.py [-h] [--work-dir WORK_DIR] --output-bids-dir OUTPUT_BIDS_DIR --stains STAINS [STAINS ...] --subject SUBJECT [--acq ACQ] [--sample SAMPLE] --input-path INPUT_PATH [--help-snakemake]

options:
  -h, --help            show this help message and exit
  --work-dir WORK_DIR   Set the work directory (effectively the snakebids output workflow directory)
  --output-bids-dir OUTPUT_BIDS_DIR
                        Set the output bids directory
  --stains STAINS [STAINS ...]
                        Set the stains for each channels
  --subject SUBJECT     Set the subject identifier (participant-label)
  --acq ACQ             Set the acquisition entity
  --sample SAMPLE       Set the sample entity
  --input-path INPUT_PATH
                        Set the input path
  --help-snakemake, --help_snakemake
                        Options to Snakemake can also be passed directly at the command-line, use this to print Snakemake usage



usage: run.py [-h] [--work-dir WORK_DIR] --output-bids-dir OUTPUT_BIDS_DIR --stains STAINS [STAINS ...] --subject SUBJECT [--acq ACQ] [--sample SAMPLE] --input-path INPUT_PATH [--help-snakemake]
Note:  the following arguments are required: --output-bids-dir, --stains, --subject, --input-path
```

e.g.: 

```
./run.py  --input-path /cifs/khan_new/datasets/lightsheet_example/ --subject M4A1Te3 --stains Abeta PI Lectin --output-bids-dir /cifs/khan_new/datasets/MIND/mouse_appmaptapoe --work-dir /tmp   --cores all --use-conda 
```


The `--input-path` should be a folder with a sample's tif files, and no other tif files, or it can point to a tar file containing the tif files.

The `--subject`, `--sample`, and `--acq` arguments identify the subject, sample, acquisition, which become part of the resulting filenames (BIDS naming).

 Note: The acquisition value must contain the strings of `blaze`, `imaris` or `prestitched`, and defines which workflow will be used. E.g. for LifeCanvas data that is already stitched, you need to include `prestitched` in the acquisition flag. 

The `--stains` argument identifies what stains are used for each channel. Standardized naming is shown in the table below (TODO: create this). 

By default the files in `--work-dir` are deleted after they are no longer needed in the workflow, unless you use the `--notemp` command-line option. The workflow writes a large number of small files in parallel to the `work` folder, so for optimum performance this should be a fast local disk, and not a networked file system (i.e. shared disk)
.

Snakemake CLI options passed to run.py are also passed along when running the workflow, to see a full list of these, use `--help-snakemake`. The required parameters from these are `--cores` and `--use-conda`, but other useful options are `--dry-run` or `-n`, and `--conda-prefix`.  

## Configuration

The `config/config.yml` can be edited to customize any workflow parameters.   



