rule extract_dataset:
    input:
        dataset_path=get_dataset_path_remote,
    params:
        cmd=cmd_extract_dataset,
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
                    suffix="SPIM",
                )
            )
        ),
    group:
        "preproc"
    log:
        bids(
            root="logs",
            subject="{subject}",
            datatype="extract_dataset",
            sample="{sample}",
            acq="{acq}",
            desc="raw",
            suffix="log.txt",
        ),
    shell:
        "{params.cmd}"


rule cp_from_gcs:
    params:
        dataset_path=get_dataset_path_gs,
    output:
        ome_dir=temp(
            directory(
                bids(
                    root=work,
                    subject="{subject}",
                    datatype="micr",
                    sample="{sample}",
                    acq="{acq}",
                    desc="rawfromgcs",
                    suffix="SPIM",
                )
            )
        ),
    threads: config["cores_per_rule"]
    group:
        "preproc"
    log:
        bids(
            root="logs",
            subject="{subject}",
            datatype="cp_from_gcs",
            sample="{sample}",
            acq="{acq}",
            desc="raw",
            suffix="log.txt",
        ),
    container:
        None
    conda:
        "../envs/google_cloud.yaml"
    shell:
        "mkdir -p {output} && gcloud storage cp --recursive {params.dataset_path}/* {output}"

rule blaze_to_metadata_gcs:
    input:
        creds = os.path.expanduser(config["remote_creds"])
    params:
        dataset_path=get_dataset_path_gs,
        in_tif_pattern=lambda wildcards: config["import_blaze"]["raw_tif_pattern"],
        storage_provider_settings=workflow.storage_provider_settings,
    output:
        metadata_json=bids(
            root=root,
            desc='gcs',
            subject="{subject}",
            datatype="micr",
            sample="{sample}",
            acq="{acq,[a-zA-Z0-9]*blaze[a-zA-Z0-9]*}",
            suffix="SPIM.json",
        ),
    benchmark:
        bids(
            root="benchmarks",
            datatype="blaze_to_metadata_gcs",
            subject="{subject}",
            sample="{sample}",
            acq="{acq}",
            suffix="benchmark.tsv",
        )
    log:
        bids(
            root="logs",
            datatype="blaze_to_metadata_gcs",
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
        "../scripts/blaze_to_metadata_gcs.py"


rule blaze_to_metadata:
    input:
        ome_dir=get_input_dataset,
    params:
        in_tif_pattern=lambda wildcards, input: os.path.join(
            input.ome_dir,
            config["import_blaze"]["raw_tif_pattern"],
        ),
    output:
        metadata_json=temp(bids(
            root=work,
            subject="{subject}",
            desc='local',
            datatype="micr",
            sample="{sample}",
            acq="{acq,[a-zA-Z0-9]*blaze[a-zA-Z0-9]*}",
            suffix="SPIM.json",
        )),
    benchmark:
        bids(
            root="benchmarks",
            datatype="blaze_to_metdata",
            subject="{subject}",
            sample="{sample}",
            acq="{acq}",
            suffix="benchmark.tsv",
        )
    log:
        bids(
            root="logs",
            datatype="blaze_to_metdata",
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
        "../scripts/blaze_to_metadata.py"


rule copy_blaze_metadata:
    input:
        json=get_metadata_json
    output:
        metadata_json=bids(
            root=root,
            subject="{subject}",
            datatype="micr",
            sample="{sample}",
            acq="{acq,[a-zA-Z0-9]*blaze[a-zA-Z0-9]*}",
            suffix="SPIM.json",
        ),
    shell: 'cp {input} {output}'

rule prestitched_to_metadata:
    input:
        ome_dir=get_input_dataset,
    params:
        physical_size_x_um=config["import_prestitched"]["physical_size_x_um"],
        physical_size_y_um=config["import_prestitched"]["physical_size_y_um"],
        physical_size_z_um=config["import_prestitched"]["physical_size_z_um"],
    output:
        metadata_json=bids(
            root=root,
            subject="{subject}",
            datatype="micr",
            sample="{sample}",
            acq="{acq,[a-zA-Z0-9]*prestitched[a-zA-Z0-9]*}",
            suffix="SPIM.json",
        ),
    benchmark:
        bids(
            root="benchmarks",
            datatype="prestitched_to_metdata",
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
        "../scripts/prestitched_to_metadata.py"


rule tif_to_zarr:
    """ use dask to load tifs in parallel and write to zarr 
        output shape is (tiles,channels,z,y,x), with the 2d 
        images as the chunks"""
    input:
        ome_dir=get_input_dataset,
        metadata_json=rules.copy_blaze_metadata.output.metadata_json,
    params:
        in_tif_pattern=lambda wildcards, input: os.path.join(
            input.ome_dir,
            config["import_blaze"]["raw_tif_pattern"],
        ),
        intensity_rescaling=config["import_blaze"]["intensity_rescaling"],
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
                    suffix="SPIM.zarr",
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
    threads: config["cores_per_rule"]
    container:
        config["containers"]["spimprep"]
    script:
        "../scripts/tif_to_zarr.py"

rule tif_to_zarr_gcs:
    """ use dask to load tifs in parallel and write to zarr 
        output shape is (tiles,channels,z,y,x), with the 2d 
        images as the chunks"""
    input:
        metadata_json=rules.copy_blaze_metadata.output.metadata_json,
        creds = os.path.expanduser(config["remote_creds"])
    params:
        dataset_path=get_dataset_path_gs,
        in_tif_pattern=lambda wildcards:
            config["import_blaze"]["raw_tif_pattern"],
        intensity_rescaling=config["import_blaze"]["intensity_rescaling"],
        storage_provider_settings=workflow.storage_provider_settings,
    output:
        zarr=temp(
            directory(
                bids(
                    root=work,
                    subject="{subject}",
                    datatype="micr",
                    sample="{sample}",
                    acq="{acq}",
                    desc="rawfromgcs",
                    suffix="SPIM.zarr",
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
    threads: config["cores_per_rule"]
    container:
        config["containers"]["spimprep"]
    script:
        "../scripts/tif_to_zarr_gcs.py"
