# Testing the workflow

## Generate testing datasets

There are two scripts that each take in a tar file of tiff files and will create a smaller version to test with. The testing scripts should be run from within the testing directory.

1. create_test_dataset

This script takes in a larger dataset and will produce a subset of it. The user can specify a given slice step and tile step in the x and y directions.

2. create_downsampled_dataset

This script can take in any size dataset and will downsample across the x, y and z. It can be run after the test dataset script is run, but the user must specify the slice step used on the first script.

## Creating the tests

Once the test datasets are created the user can then generate the unit tests for the workflow, following these steps:

1. Change the following relative paths to absolute path.

	1. The path to config file from within the snakefile
	
	2. The path to the datasets.tsv file from within the config file
	
	3. The path to the test dataset from within the datasets.tsv file
	
	Making thesse changes will ensure the unit test have all relevant context
	
2. Run the generate_test python script from the spimprep directory with:

  ```
  python testing/generate_test.py
  ```
  
  This will run the snakemake workflow, generate the unit test and then copy in the modified test scripts to make sure the tests are correct.
  
