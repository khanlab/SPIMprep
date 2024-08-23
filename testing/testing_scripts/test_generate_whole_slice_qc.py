import os
import sys

import subprocess as sp
from tempfile import TemporaryDirectory
import shutil
from pathlib import Path, PurePosixPath

sys.path.insert(0, os.path.dirname(__file__))

import common


def test_generate_whole_slice_qc():

    with TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir) / "workdir"
        data_path = PurePosixPath(".tests/unit/generate_whole_slice_qc/data")
        expected_path = PurePosixPath(".tests/unit/generate_whole_slice_qc/expected")
        qc = PurePosixPath("qc/resources")

        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir)
        shutil.copytree(qc, workdir / "qc" / "resources")

        # dbg
        print("qc/sub-mouse1_sample-brain_acq-blaze1x/whole_slice_qc.html qc/sub-mouse1_sample-brain_acq-blaze1x/images/whole", file=sys.stderr)

        # Run the test job.
        sp.check_output([
            "python",
            "-m",
            "snakemake", 
            "qc/sub-mouse1_sample-brain_acq-blaze1x/whole_slice_qc.html",
            "qc/sub-mouse1_sample-brain_acq-blaze1x/images/whole",
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
        common.OutputChecker(data_path, expected_path, workdir).check()
