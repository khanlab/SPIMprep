# Workflow

Under construction!

## Overview

SPIMprep is implemented in Python using the [Snakemake](https://snakemake.github.io/) workflow management system. It performs metadata extraction, flatfield correction (BaSiC), and stitching (BigStitcher), followed by the creation of a final validated BIDS dataset. 

## Inputs to SPIMprep




## Installation

The workflow is installed via pip, and the required container dependency is downloaded automatically by Snakemake. 

## Configuration

Input datasets are configured using a TSV file, specifying subject identifiers and paths to folders/archives containing the raw TIF files, and a YAML file is used to customize workflow configuration. 

## Running Parallelization

The workflow can be executed in parallel on any local, cluster or cloud resources, and each step is also internally parallelized (Dask), taking advantage of the parallelization afforded by the chunked file format. 

## Cloud support

The workflow can optionally write directly to cloud storage, facilitating data sharing and interoperability with existing web-based viewers. 

## Outputs of SPIMprep


