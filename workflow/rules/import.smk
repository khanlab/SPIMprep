
rule get_dataset:
    input:
        dataset_path=get_dataset_path,
    params:
        cmd=cmd_get_dataset,
    output:
        ome_dir=temp(
            directory(
                bids(
                    root=work,
                    subject="{subject}",
                    datatype="micr",
                    sample="{sample}",
                    acq="{acq}",
                    desc="raw",
                    suffix="spim",
                )
            )
        ),
    group:
        "preproc"
    log:
        bids(
            root="logs",
            subject="{subject}",
            datatype="get_dataset",
            sample="{sample}",
            acq="{acq}",
            desc="raw",
            suffix="log.txt",
        ),
    shell:
        "{params.cmd}"


rule raw_to_metadata:
    input:
        ome_dir=rules.get_dataset.output.ome_dir,
    params:
        in_tif_glob=lambda wildcards, input: os.path.join(
            input.ome_dir,
            config["import"]["raw_tif_pattern"],
        ),
    output:
        metadata_json=bids(
            root=root,
            subject="{subject}",
            datatype="micr",
            sample="{sample}",
            acq="{acq}",
            suffix="spim.metadata.json",
        ),
    benchmark:
        bids(
            root="benchmarks",
            datatype="raw_to_metdata",
            subject="{subject}",
            sample="{sample}",
            acq="{acq}",
            suffix="benchmark.tsv",
        )
    log:
        bids(
            root="logs",
            datatype="raw_to_metdata",
            subject="{subject}",
            sample="{sample}",
            acq="{acq}",
            suffix="log.txt",
        ),
    group:
        "preproc"
    container:
        config["containers"]["spimprep"]
    script:
        "../scripts/raw_to_metadata.py"


rule tif_to_zarr:
    """ use dask to load tifs in parallel and write to zarr 
        output shape is (tiles,channels,z,y,x), with the 2d 
        images as the chunks"""
    input:
        ome_dir=rules.get_dataset.output.ome_dir,
        metadata_json=rules.raw_to_metadata.output.metadata_json,
    params:
        in_tif_glob=lambda wildcards, input: os.path.join(
            input.ome_dir,
            config["import"]["raw_tif_glob"],
        ),
        intensity_rescaling=config["import"]["intensity_rescaling"],
    output:
        zarr=temp(
            directory(
                bids(
                    root=work,
                    subject="{subject}",
                    datatype="micr",
                    sample="{sample}",
                    acq="{acq}",
                    desc="raw",
                    suffix="spim.zarr",
                )
            )
        ),
    benchmark:
        bids(
            root="benchmarks",
            datatype="tif_to_zarr",
            subject="{subject}",
            sample="{sample}",
            acq="{acq}",
            suffix="benchmark.tsv",
        )
    log:
        bids(
            root="logs",
            datatype="tif_to_zarr",
            subject="{subject}",
            sample="{sample}",
            acq="{acq}",
            suffix="log.txt",
        ),
    group:
        "preproc"
    threads: 32
    container:
        config["containers"]["spimprep"]
    script:
        "../scripts/tif_to_zarr.py"
