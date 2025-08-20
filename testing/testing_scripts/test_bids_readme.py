import os
import sys

import subprocess as sp
from tempfile import TemporaryDirectory
import shutil
from pathlib import Path, PurePosixPath

sys.path.insert(0, os.path.dirname(__file__))

import common


def test_bids_readme():

    with TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir) / "workdir"
        data_path = PurePosixPath(".tests/unit/bids_readme/data")
        expected_path = PurePosixPath(".tests/unit/bids_readme/expected")

        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir)

        # dbg
        print("bids/README.md", file=sys.stderr)
        #Hello

        # Run the test job.
        sp.check_output([
            "python",
            "-m",
            "snakemake", 
            "bids/README.md",
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
