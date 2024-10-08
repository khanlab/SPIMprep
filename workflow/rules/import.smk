rule extract_sample:
    input:
        sample_path=get_sample_path_remote,
    params:
        cmd=cmd_extract_sample,
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
            datatype="extract_sample",
            sample="{sample}",
            acq="{acq}",
            desc="raw",
            suffix="log.txt",
        ),
    shell:
        "{params.cmd}"


rule blaze_to_metadata_gcs:
    input:
        creds=os.path.expanduser(config["remote_creds"]),
    params:
        sample_path=get_sample_path_gs,
        in_tif_pattern=lambda wildcards: config["import_blaze"]["raw_tif_pattern"],
        storage_provider_settings=workflow.storage_provider_settings,
    output:
        metadata_json=bids(
            root=work,
            desc="gcs",
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
        ome_dir=get_input_sample,
    output:
        metadata_json=temp(
            bids(
                root=work,
                subject="{subject}",
                desc="local",
                datatype="micr",
                sample="{sample}",
                acq="{acq,[a-zA-Z0-9]*blaze[a-zA-Z0-9]*}",
                suffix="SPIM.json",
            )
        ),
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
        json=get_metadata_json,
    output:
        metadata_json=bids(
            root=root,
            subject="{subject}",
            datatype="micr",
            sample="{sample}",
            acq="{acq,[a-zA-Z0-9]*blaze[a-zA-Z0-9]*}",
            suffix="SPIM.json",
        ),
    log:
        bids(
            root="logs",
            datatype="copy_blaze_metadata",
            subject="{subject}",
            sample="{sample}",
            acq="{acq}",
            suffix="log.txt",
        ),
    shell:
        "cp {input} {output} &> {log}"


rule prestitched_to_metadata:
    input:
        ome_dir=get_input_sample,
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
        ome_dir=get_input_sample,
        metadata_json=rules.copy_blaze_metadata.output.metadata_json,
    params:
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
    resources:
        mem_mb=config["total_mem_mb"],
    threads: int(config["total_mem_mb"] / 8000)  #this is memory-limited -- seems to need ~8000mb for each thread, so threads=total_mem_mb / 8000 
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
        creds=os.path.expanduser(config["remote_creds"]),
    params:
        sample_path=get_sample_path_gs,
        in_tif_pattern=lambda wildcards: config["import_blaze"]["raw_tif_pattern"],
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
    resources:
        mem_mb=config["total_mem_mb"],
    threads: int(config["total_mem_mb"] / 8000)  #this is memory-limited -- seems to need ~8000mb for each thread, so threads=total_mem_mb / 8000 
    container:
        config["containers"]["spimprep"]
    script:
        "../scripts/tif_to_zarr_gcs.py"
