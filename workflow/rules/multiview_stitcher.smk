ruleorder: multiview_stitcher > zarr_to_ome_zarr

rule multiview_stitcher:
    input:
        **get_storage_creds(),
        zarr=bids(
            root=work,
            subject="{subject}",
            datatype="micr",
            sample="{sample}",
            acq="{acq}",
            desc=config["ome_zarr"]["desc"],
            suffix="SPIM.zarr",
        ),
        metadata_json=rules.copy_blaze_metadata.output.metadata_json,
    params:
        channels=get_stains,
        registration_opts=config['multiview_stitcher']['registration'],
        fusion_opts=config['multiview_stitcher']['fusion'],
        reg_channel_index=1, #later can make this a parameter chosen based on stains
        uri=get_output_ome_zarr_uri(),
        storage_provider_settings=workflow.storage_provider_settings,
    output:
        **get_output_ome_zarr("blaze"),
#    shadow: 'minimal'  #don't make this shadow, as we then have to copy files after (less efficient)
    #instead just make a work temp() folder for the intermediate zarrs - better for cloud too..
    benchmark:
        bids(
            root="benchmarks",
            datatype="multiview_stitcher",
            subject="{subject}",
            sample="{sample}",
            acq="{acq}",
            suffix="benchmark.tsv",
        )
    log:
        bids(
            root="logs",
            datatype="multiview_stitcher",
            subject="{subject}",
            sample="{sample}",
            acq="{acq}",
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


