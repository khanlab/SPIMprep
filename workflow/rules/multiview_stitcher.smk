ruleorder: multiview_stitcher > zarr_to_ome_zarr

rule multiview_stitcher:
    #note: this runs once for each channel (though the registration should be identical
    input:
        **get_storage_creds(),
        zarr=bids(
            root=work,
            subject="{subject}",
            datatype="micr",
            sample="{sample}",
            acq="{acq}",
            desc="{desc}",
            suffix="SPIM.zarr",
        ),
        metadata_json=rules.copy_blaze_metadata.output.metadata_json,
    params:
        channels=get_stains,
        channel_index=lambda wildcards: get_stains(wildcards).index(wildcards.stain),
        registration_opts=config['multiview_stitcher']['registration'],
        fusion_opts=config['multiview_stitcher']['fusion'],
        reg_channel_index=1, #later can make this a parameter chosen based on stains
        uri=get_output_ome_zarr_uri(),
        storage_provider_settings=workflow.storage_provider_settings,
    output:
        zarr=directory(bids(
            root=work,
            subject="{subject}",
            datatype="micr",
            sample="{sample}",
            acq="{acq}",
            desc="mvstitched{desc}",
            stain="{stain}",
            suffix="SPIM.zarr",
        )),
    benchmark:
        bids(
            root="benchmarks",
            datatype="multiview_stitcher",
            subject="{subject}",
            sample="{sample}",
            acq="{acq}",
            desc="{desc}",
            stain="{stain}",
            suffix="benchmark.tsv",
        )
    log:
        bids(
            root="logs",
            datatype="multiview_stitcher",
            subject="{subject}",
            sample="{sample}",
            acq="{acq}",
            desc="{desc}",
            stain="{stain}",
            suffix="log.txt",
        ),
    threads: config["total_cores"]
    resources:
        mem_mb=config["total_mem_mb"],
    group:
        "preproc"
    container: None
    conda:
        "../envs/multiview_stitcher.yml"
    script:
        "../scripts/multiview_stitcher.py"


