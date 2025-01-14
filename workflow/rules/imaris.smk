rule imaris_to_metadata:
    input:
        ims=get_input_sample,
    output:
        metadata_json=bids(
            root=root,
            subject="{subject}",
            datatype="micr",
            sample="{sample}",
            acq="{acq,[a-zA-Z0-9]*imaris[a-zA-Z0-9]*}",
            suffix="SPIM.json",
        ),
    benchmark:
        bids(
            root="benchmarks",
            datatype="imaris_to_metdata",
            subject="{subject}",
            sample="{sample}",
            acq="{acq}",
            suffix="benchmark.tsv",
        )
    log:
        bids(
            root="logs",
            datatype="prestitched_to_metdata",
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
        "../scripts/imaris_to_metadata.py"

rule imaris_channel_to_zarr:
    input:
        ims=get_input_sample,
    params:
        channel=lambda wildcards: get_stains(wildcards).index(wildcards.stain),
    output:
        zarr=directory(bids(
            root=root,
            subject="{subject}",
            datatype="micr",
            sample="{sample}",
            acq="{acq}",
            stain="{stain}",
            suffix="imaris.zarr",
        )),
    log:
        bids(
            root="logs",
            subject="{subject}",
            datatype="imaris_channel_to_zarr",
            sample="{sample}",
            acq="{acq}",
            stain="{stain}",
            suffix="log.txt",
        ),
    container:
        config["containers"]["spimprep"]
    group:
        "preproc"
    threads: 1
    resources:
        runtime=360,
        mem_mb=1000,
    shadow: 'minimal'
    script:
        "../scripts/imaris_channel_to_zarr.py"


rule imaris_to_ome_zarr:
    input:
        ims=get_input_sample,
        zarr=lambda wildcards: expand(bids(
            root=root,
            subject="{subject}",
            datatype="micr",
            sample="{sample}",
            acq="{acq}",
            stain="{stain}",
            suffix="imaris.zarr",
        ),stain=get_stains(wildcards),allow_missing=True),
        metadata_json=rules.prestitched_to_metadata.output.metadata_json,
    params:
        max_downsampling_layers=config["ome_zarr"]["max_downsampling_layers"],
        rechunk_size=config["ome_zarr"]["rechunk_size"],
        scaling_method=config["ome_zarr"]["scaling_method"],
        downsampling=config["bigstitcher"]["fuse_dataset"]["downsampling"],
        stains=get_stains,
        uri=get_output_ome_zarr_uri(),
        storage_provider_settings=workflow.storage_provider_settings,
    output:
        **get_output_ome_zarr("imaris"),
    log:
        bids(
            root="logs",
            subject="{subject}",
            datatype="imaris_to_ome_zarr",
            sample="{sample}",
            acq="{acq}",
            suffix="log.txt",
        ),
    container:
        config["containers"]["spimprep"]
    group:
        "preproc"
    threads: config["total_cores"]
    resources:
        runtime=360,
        mem_mb=config["total_mem_mb"],
    shadow: 'minimal'
    script:
        "../scripts/imaris_to_ome_zarr.py"


