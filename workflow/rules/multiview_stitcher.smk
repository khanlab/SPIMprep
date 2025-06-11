
rule mvstitcher_registration:
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
        registration_opts=config['multiview_stitcher']['registration'],
        fusion_opts=config['multiview_stitcher']['fusion'],
        reg_channel_index=1, #later can make this a parameter chosen based on stains
        uri=get_output_ome_zarr_uri(),
        storage_provider_settings=workflow.storage_provider_settings,
    output:
        affines=bids(
            root=work,
            subject="{subject}",
            datatype="micr",
            sample="{sample}",
            acq="{acq}",
            desc="mvstitched{desc}",
            suffix="affines.npz",
        ),
    benchmark:
        bids(
            root="benchmarks",
            datatype="mvstitcher_registration",
            subject="{subject}",
            sample="{sample}",
            acq="{acq}",
            desc="{desc}",
            suffix="benchmark.tsv",
        )
    log:
        bids(
            root="logs",
            datatype="mvstitcher_registration",
            subject="{subject}",
            sample="{sample}",
            acq="{acq}",
            desc="{desc}",
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
        "../scripts/mvstitcher_registration.py"

ruleorder: mvstitcher_fusion > zarr_to_ome_zarr

rule mvstitcher_fusion:
    input:
        **get_storage_creds(),
        zarr=bids(
            root=work,
            subject="{subject}",
            datatype="micr",
            sample="{sample}",
            acq="{acq}",
            desc="flatcorr",
            suffix="SPIM.zarr",
        ),
        metadata_json=rules.copy_blaze_metadata.output.metadata_json,
        affines=bids(
            root=work,
            subject="{subject}",
            datatype="micr",
            sample="{sample}",
            acq="{acq}",
            desc="mvstitchedflatcorr",
            suffix="affines.npz",
        ),
    params:
        channels=get_stains,
        registration_opts=config['multiview_stitcher']['registration'],
        fusion_opts=config['multiview_stitcher']['fusion'],
        reg_channel_index=1, #later can make this a parameter chosen based on stains
        uri=get_output_ome_zarr_uri(),
        storage_provider_settings=workflow.storage_provider_settings,
    output:
        **get_output_ome_zarr("blaze"),
    benchmark:
        bids(
            root="benchmarks",
            datatype="mvstitcher_fusion",
            subject="{subject}",
            sample="{sample}",
            acq="{acq}",
            suffix="benchmark.tsv",
        )
    log:
        bids(
            root="logs",
            datatype="mvstitcher_fusion",
            subject="{subject}",
            sample="{sample}",
            acq="{acq}",
            suffix="log.txt",
        ),
    threads: 128
    resources:
        mem_mb=config["total_mem_mb"],
    group:
        "preproc"
    container: None
    conda:
        "../envs/multiview_stitcher.yml"
    script:
        "../scripts/mvstitcher_fusion.py"


