

rule zarr_to_ome_zarr:
    input:
        **get_storage_creds(),
        zarr=lambda wildcards: expand(
            bids(
                root=work,
                subject="{subject}",
                datatype="micr",
                sample="{sample}",
                acq="{acq}",
                desc="{desc}",
                stain="{stain}",
                suffix="SPIM.zarr",
            ),
            stain=get_stains(wildcards),
            desc=config["ome_zarr"]["desc"],
            allow_missing=True,
        ),
        metadata_json=rules.copy_blaze_metadata.output.metadata_json,
    params:
        max_downsampling_layers=config["ome_zarr"]["max_downsampling_layers"],
        rechunk_size=config["ome_zarr"]["rechunk_size"],
        scaling_method=config["ome_zarr"]["scaling_method"],
        downsampling=config["bigstitcher"]["fuse_dataset"]["downsampling"],
        stains=get_stains,
        uri=get_output_ome_zarr_uri(),
        storage_provider_settings=workflow.storage_provider_settings,
    output:
        **get_output_ome_zarr("blaze"),
    resources:
        mem_mb=config["total_mem_mb"],
    threads: config["total_cores"]
    log:
        bids(
            root="logs",
            subject="{subject}",
            datatype="zarr_to_ome_zarr",
            sample="{sample}",
            acq="{acq}",
            suffix="log.txt",
        ),
    container:
        config["containers"]["spimprep"]
    group:
        "preproc"
    script:
        "../scripts/zarr_to_ome_zarr.py"


rule tif_stacks_to_ome_zarr:
    input:
        **get_storage_creds(),
        tif_dir=get_input_sample,
        metadata_json=rules.prestitched_to_metadata.output.metadata_json,
    params:
        in_tif_glob=lambda wildcards, input: os.path.join(
            input.tif_dir,
            config["import_prestitched"]["stitched_tif_glob"],
        ),
        max_downsampling_layers=config["ome_zarr"]["max_downsampling_layers"],
        rechunk_size=config["ome_zarr"]["rechunk_size"],
        scaling_method=config["ome_zarr"]["scaling_method"],
        downsampling=config["bigstitcher"]["fuse_dataset"]["downsampling"],
        stains=get_stains,
        uri=get_output_ome_zarr_uri(),
        storage_provider_settings=workflow.storage_provider_settings,
    output:
        **get_output_ome_zarr("prestitched"),
    log:
        bids(
            root="logs",
            subject="{subject}",
            datatype="tif_stacks_to_ome_zarr",
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
    script:
        "../scripts/tif_stacks_to_ome_zarr.py"



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
    container:
        config["containers"]["spimprep"]
    script:
        "../scripts/ome_zarr_to_nii.py"

