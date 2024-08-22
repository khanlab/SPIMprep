import os
import sys

import subprocess as sp
from tempfile import TemporaryDirectory
import shutil
from pathlib import Path, PurePosixPath

sys.path.insert(0, os.path.dirname(__file__))

import common


def test_generate_flatfield_qc():

    with TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir) / "workdir"
        data_path = PurePosixPath(".tests/unit/generate_flatfield_qc/data")
        expected_path = PurePosixPath(".tests/unit/generate_flatfield_qc/expected")
        qc = PurePosixPath("qc/resources")

        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir)
        shutil.copytree(qc, workdir / "qc" / "resources")

        # dbg
        print("qc/sub-mouse1_sample-brain_acq-blaze1x/flatfieldqc.html qc/sub-mouse1_sample-brain_acq-blaze1x/images/corr qc/sub-mouse1_sample-brain_acq-blaze1x/images/uncorr", file=sys.stderr)

        # Run the test job.
        sp.check_output([
            "python",
            "-m",
            "snakemake", 
            "qc/sub-mouse1_sample-brain_acq-blaze1x/flatfieldqc.html",
            "qc/sub-mouse1_sample-brain_acq-blaze1x/images/corr",
            "qc/sub-mouse1_sample-brain_acq-blaze1x/images/uncorr",
            "-f", 
            "-j1",
            "--keep-target-files",
    
            "--directory",
            workdir,
        ])

        # Check the output byte by byte using cmp.
        # To modify this behavior, you can inherit from common.OutputChecker in here
        # and overwrite the method `compare_files(generated_file, expected_file), 
        # also see common.py.
        common.OutputChecker(data_path, expected_path, workdir).check()
