import subprocess as sp
import os
from pathlib import Path
import shutil

# Run snakemake workflow with no temp to be able to generate tests
sp.run([
    "python",
    "-m",
    "snakemake", 
    "-c",
    "all", 
    "--use-singularity",
    "--notemp"
])

# Generate the unit tests
sp.run([
    "python",
    "-m",
    "snakemake",
    "--generate-unit-tests"
])

# Path to testing scripts
directory = Path("testing/testing_scripts")
# output the tests in the unit test folder
output_directory = Path(".tests/unit")
# get all the test files
files = os.listdir(directory)

# Copy all the test scripts into the unit test directory
for file in files:
    full_name = directory / file
    full_output_name = output_directory / file
    shutil.copy(full_name, full_output_name)




