import os
import sys

import subprocess as sp
from tempfile import TemporaryDirectory
import shutil
from pathlib import Path, PurePosixPath

sys.path.insert(0, os.path.dirname(__file__))

import common


def test_zarr_to_bdv():

    with TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir) / "workdir"
        data_path = PurePosixPath(".tests/unit/zarr_to_bdv/data")
        expected_path = PurePosixPath(".tests/unit/zarr_to_bdv/expected")

        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir)

        # dbg
        print("work/sub-mouse1/micr/sub-mouse1_sample-brain_acq-blaze1x_desc-flatcorr_bdv.n5 work/sub-mouse1/micr/sub-mouse1_sample-brain_acq-blaze1x_desc-flatcorr_bdv.xml", file=sys.stderr)

        # Run the test job.
        sp.check_output([
            "python",
            "-m",
            "snakemake", 
            "work/sub-mouse1/micr/sub-mouse1_sample-brain_acq-blaze1x_desc-flatcorr_bdv.n5",
            "work/sub-mouse1/micr/sub-mouse1_sample-brain_acq-blaze1x_desc-flatcorr_bdv.xml",
            "-f", 
            "-j1",
            "--target-files-omit-workdir-adjustment",
			"--use-singularity",
    
            "--directory",
            workdir,
        ])

        # Check the output byte by byte using cmp.
        # To modify this behavior, you can inherit from common.OutputChecker in here
        # and overwrite the method `compare_files(generated_file, expected_file), 
        # also see common.py.
        common.ImperfectOutputChecker(data_path, expected_path, workdir).check()
