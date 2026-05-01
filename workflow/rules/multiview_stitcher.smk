
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
        reg_result_pkl=bids(
            root=work,
            subject="{subject}",
            datatype="micr",
            sample="{sample}",
            acq="{acq}",
            desc="mvstitched{desc}",
            suffix="regresult.pkl",
        ),
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
    conda:
        "../envs/multiview_stitcher.yml"
    container:
        None
    threads: config["total_cores"]
    resources:
        mem_mb=config["total_mem_mb"],
        runtime=240,
    params:
        channels=get_stains,
        registration_opts=config["multiview_stitcher"]["registration"],
        fusion_opts=config["multiview_stitcher"]["fusion"],
        reg_channel_index=1,  #later can make this a parameter chosen based on stains
        uri=get_output_ome_zarr_uri(),
        storage_provider_settings=workflow.storage_provider_settings,
    script:
        "../scripts/mvstitcher_registration.py"


rule mvstitcher_reg_plots:
    input:
        reg_result_pkl=bids(
            root=work,
            subject="{subject}",
            datatype="micr",
            sample="{sample}",
            acq="{acq}",
            desc="mvstitched{desc}",
            suffix="regresult.pkl",
        ),
    output:
        pairwise_plot=bids(
            root=work,
            subject="{subject}",
            datatype="micr",
            sample="{sample}",
            acq="{acq}",
            desc="mvstitched{desc}",
            suffix="pairwiseqc.png",
        ),
        groupwise_plot=bids(
            root=work,
            subject="{subject}",
            datatype="micr",
            sample="{sample}",
            acq="{acq}",
            desc="mvstitched{desc}",
            suffix="groupwiseqc.png",
        ),
    conda:
        "../envs/multiview_stitcher.yml"
    container:
        None
    threads: 1
    resources:
        mem_mb=4000,
        runtime=30,
    script:
        "../scripts/mvstitcher_reg_plots.py"


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
            desc="{desc}".format(
                desc="flatcorr" if config["basic_flatfield_corr"]["enabled"] else "raw"
            ),
            suffix="SPIM.zarr",
        ),
        metadata_json=rules.copy_blaze_metadata.output.metadata_json,
        affines=bids(
            root=work,
            subject="{subject}",
            datatype="micr",
            sample="{sample}",
            acq="{acq}",
            desc="mvstitched{desc}".format(
                desc="flatcorr" if config["basic_flatfield_corr"]["enabled"] else "raw"
            ),
            suffix="affines.npz",
        ),
    output:
        **get_output_ome_zarr("blaze"),
    log:
        bids(
            root="logs",
            datatype="mvstitcher_fusion",
            subject="{subject}",
            sample="{sample}",
            acq="{acq}",
            suffix="log.txt",
        ),
    benchmark:
        bids(
            root="benchmarks",
            datatype="mvstitcher_fusion",
            subject="{subject}",
            sample="{sample}",
            acq="{acq}",
            suffix="benchmark.tsv",
        )
    conda:
        "../envs/multiview_stitcher.yml"
    container:
        None
    threads: 128
    resources:
        mem_mb=config["total_mem_mb"],
        runtime=480,
    params:
        channels=get_stains,
        registration_opts=config["multiview_stitcher"]["registration"],
        fusion_opts=config["multiview_stitcher"]["fusion"],
        reg_channel_index=1,  #later can make this a parameter chosen based on stains
        uri=get_output_ome_zarr_uri(),
        storage_provider_settings=workflow.storage_provider_settings,
    script:
        "../scripts/mvstitcher_fusion.py"
