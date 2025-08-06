rule ome_zarr_to_nii:
    input:
        **get_storage_creds(),
        zarr=get_input_ome_zarr_to_nii,
    params:
        channel_index=lambda wildcards: get_stains(wildcards).index(wildcards.stain),
        uri=get_output_ome_zarr_uri(),
        storage_provider_settings=workflow.storage_provider_settings,
    output:
        nii=bids(
            root=resampled,
            subject="{subject}",
            datatype="micr",
            sample="{sample}",
            acq="{acq}",
            res="{level}x",
            stain="{stain}",
            suffix="SPIM.nii",
        ),
    benchmark:
        bids(
            root="benchmarks",
            datatype="ome_zarr_to_nifti",
            subject="{subject}",
            sample="{sample}",
            acq="{acq}",
            res="{level}x",
            stain="{stain}",
            suffix="benchmark.tsv",
        )
    log:
        bids(
            root="logs",
            datatype="ome_zarr_to_nifti",
            subject="{subject}",
            sample="{sample}",
            acq="{acq}",
            res="{level}x",
            stain="{stain}",
            suffix="log.txt",
        ),
    group:
        "preproc"
    threads: config["total_cores"]
    container: None
#        config["containers"]["spimprep"]
    script:
        "../scripts/ome_zarr_to_nii.py"
