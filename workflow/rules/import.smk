rule extract_sample:
    input:
        sample_path=get_sample_path_remote,
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
    group:
        "preproc"
    threads: 1
    resources:
        mem_mb=4000,
        runtime=60,
    params:
        cmd=cmd_extract_sample,
    shell:
        "{params.cmd}"


rule blaze_to_metadata_gcs:
    input:
        creds=os.path.expanduser(config["remote_creds"]),
    output:
        metadata_json=bids(
            root=work,
            desc="gcs",
            subject="{subject}",
            datatype="micr",
            sample="{sample}",
            acq="{acq}",
            suffix="SPIM.json",
        ),
    log:
        bids(
            root="logs",
            datatype="blaze_to_metadata_gcs",
            subject="{subject}",
            sample="{sample}",
            acq="{acq}",
            suffix="log.txt",
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
    wildcard_constraints:
        acq=get_acq_constraint("blaze"),
    group:
        "preproc"
    container:
        config["containers"]["spimprep"]
    threads: 1
    resources:
        mem_mb=2000,
        runtime=60,
    params:
        sample_path=get_sample_path_gs,
        in_tif_pattern=lambda wildcards: config["import_blaze"]["raw_tif_pattern"],
        storage_provider_settings=workflow.storage_provider_settings,
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
                acq="{acq}",
                suffix="SPIM.json",
            )
        ),
    log:
        bids(
            root="logs",
            datatype="blaze_to_metdata",
            subject="{subject}",
            sample="{sample}",
            acq="{acq}",
            suffix="log.txt",
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
    wildcard_constraints:
        acq=get_acq_constraint("blaze"),
    group:
        "preproc"
    container:
        config["containers"]["spimprep"]
    threads: 1
    resources:
        mem_mb=2000,
        runtime=60,
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
            acq="{acq}",
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
    wildcard_constraints:
        acq=get_acq_constraint("blaze"),
    threads: 1
    resources:
        mem_mb=1000,
        runtime=10,
    shell:
        "cp {input} {output} &> {log}"


rule prestitched_to_metadata:
    input:
        ome_dir=get_input_sample,
    output:
        metadata_json=bids(
            root=root,
            subject="{subject}",
            datatype="micr",
            sample="{sample}",
            acq="{acq}",
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
            datatype="prestitched_to_metdata",
            subject="{subject}",
            sample="{sample}",
            acq="{acq}",
            suffix="benchmark.tsv",
        )
    wildcard_constraints:
        acq=get_acq_constraint("prestitched"),
    group:
        "preproc"
    container:
        config["containers"]["spimprep"]
    threads: 1
    resources:
        mem_mb=2000,
        runtime=60,
    params:
        physical_size_x_um=config["import_prestitched"]["physical_size_x_um"],
        physical_size_y_um=config["import_prestitched"]["physical_size_y_um"],
        physical_size_z_um=config["import_prestitched"]["physical_size_z_um"],
    script:
        "../scripts/prestitched_to_metadata.py"


rule tif_to_zarr_gcs:
    """use dask to load tifs in parallel and write to zarr
output shape is (tiles,channels,z,y,x), with the 2d
images as the chunks"""
    input:
        metadata_json=rules.copy_blaze_metadata.output.metadata_json,
        creds=os.path.expanduser(config["remote_creds"]),
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
    log:
        bids(
            root="logs",
            datatype="tif_to_zarr",
            subject="{subject}",
            sample="{sample}",
            acq="{acq}",
            suffix="log.txt",
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
    group:
        "preproc"
    container:
        config["containers"]["spimprep"]
    threads: int(config["total_mem_mb"] / 8000)  # this is memory-limited -- seems to need ~8000mb for each thread, so threads=total_mem_mb / 8000 
    resources:
        mem_mb=config["total_mem_mb"],
        runtime=360,
    params:
        sample_path=get_sample_path_gs,
        in_tif_pattern=lambda wildcards: config["import_blaze"]["raw_tif_pattern"],
        intensity_rescaling=config["import_blaze"]["intensity_rescaling"],
        storage_provider_settings=workflow.storage_provider_settings,
    script:
        "../scripts/tif_to_zarr_gcs.py"


rule bioformats_to_zarr:
    """
Use bioformats2raw on each tile, then put all tiles into a single zarr dataset.
Output shape is (tiles,channels,z,y,x), with the 2D images as the chunks.
TODO: this could potentially be done in parallel, e.g. using wildcards over tile identifiers.

"""
    input:
        [],
    output:
        tiles_dir=temp(
            directory(
                bids(
                    root=work,
                    subject="{subject}",
                    datatype="micr",
                    sample="{sample}",
                    acq="{acq}",
                    desc="raw",
                    suffix="SPIM.tiles",
                )
            )
        ),
    log:
        bids(
            root="logs",
            datatype="bioformats_to_zarr",
            subject="{subject}",
            sample="{sample}",
            acq="{acq}",
            suffix="log.txt",
        ),
    benchmark:
        bids(
            root="benchmarks",
            datatype="bioformats_to_zarr",
            subject="{subject}",
            sample="{sample}",
            acq="{acq}",
            suffix="benchmark.tsv",
        )
    group:
        "preproc"
    threads: 16
    resources:
        mem_mb=config["total_mem_mb"],  # TODO update this, along with threads.. 
        runtime=360,
        disk_mb=1000000,  #1TB
    #       ome_dir=get_input_sample
    params:
        ome_dir=get_input_sample,
        tile_height=4096,
        tile_width=4096,
    script:
        "../scripts/bioformats_to_zarr.py"


rule concat_tiles:
    """read in zarrs created for each tile, and write out as a single zarr"""
    input:
        tiles_dir=bids(
            root=work,
            subject="{subject}",
            datatype="micr",
            sample="{sample}",
            acq="{acq}",
            desc="raw",
            suffix="SPIM.tiles",
        ),
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
    log:
        bids(
            root="logs",
            datatype="concat_tiles",
            subject="{subject}",
            sample="{sample}",
            acq="{acq}",
            suffix="log.txt",
        ),
    benchmark:
        bids(
            root="benchmarks",
            datatype="concat_tiles",
            subject="{subject}",
            sample="{sample}",
            acq="{acq}",
            suffix="benchmark.tsv",
        )
    group:
        "preproc"
    container:
        None
    threads: 32
    resources:
        mem_mb=config["total_mem_mb"],
        runtime=240,
        disk_mb=1000000,  #1TB
    params:
        intensity_rescaling=config["import_blaze"]["intensity_rescaling"],
    script:
        "../scripts/concat_tiles.py"
