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
    log:
        bids(
            root="logs",
            datatype="prestitched_to_metdata",
            subject="{subject}",
            sample="{sample}",
            acq="{acq}",
            suffix="log.txt",
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
    conda:
        "../envs/imaris.yml"
    threads: 1
    resources:
        mem_mb=2000,
        runtime=60,
    script:
        "../scripts/imaris_to_metadata.py"


rule imaris_channel_to_zarr:
    input:
        ims=get_input_sample,
    output:
        zarr=temp(
            directory(
                bids(
                    root=work,
                    subject="{subject}",
                    datatype="micr",
                    sample="{sample}",
                    acq="{acq}",
                    stain="{stain}",
                    suffix="imaris.zarr",
                )
            ),
            group_jobs=True,
        ),
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
    shadow:
        "minimal"
    threads: 1
    resources:
        runtime=360,
        mem_mb=32000,
    params:
        channel=lambda wildcards: get_stains(wildcards).index(wildcards.stain),
    script:
        "../scripts/imaris_channel_to_zarr.py"


rule imaris_to_ome_zarr:
    input:
        zarr=lambda wildcards: expand(
            bids(
                root=work,
                subject="{subject}",
                datatype="micr",
                sample="{sample}",
                acq="{acq}",
                stain="{stain}",
                suffix="imaris.zarr",
            ),
            stain=get_stains(wildcards),
            allow_missing=True,
        ),
        metadata_json=rules.prestitched_to_metadata.output.metadata_json,
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
    shadow:
        "minimal"
    threads: config["total_cores"]
    resources:
        runtime=360,
        mem_mb=300000,
    params:
        max_downsampling_layers=config["ome_zarr"]["max_downsampling_layers"],
        rechunk_size=config["ome_zarr"]["rechunk_size"],
        scaling_method=config["ome_zarr"]["scaling_method"],
        downsampling=config["ome_zarr"]["downsampling"],
        stains=get_stains,
        uri=get_output_ome_zarr_uri(),
        storage_provider_settings=workflow.storage_provider_settings,
    script:
        "../scripts/imaris_to_ome_zarr.py"
