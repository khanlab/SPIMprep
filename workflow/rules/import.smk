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


rule bioformats_to_zarr:
    """
    Use bioformats2raw on each tile, then put all tiles into a single zarr dataset.
    Output shape is (tiles,channels,z,y,x), with the 2D images as the chunks.
    TODO: this could potentially be done in parallel, e.g. using wildcards over tile identifiers.
    """
    input:
        "ome_dir=get_input_sample",
    #        ome_dir=get_input_sample
    params:
        ome_dir=get_input_sample,
        tile_height=4096,
        tile_width=4096,
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
    benchmark:
        bids(
            root="benchmarks",
            datatype="bioformats_to_zarr",
            subject="{subject}",
            sample="{sample}",
            acq="{acq}",
            suffix="benchmark.tsv",
        )
    log:
        bids(
            root="logs",
            datatype="bioformats_to_zarr",
            subject="{subject}",
            sample="{sample}",
            acq="{acq}",
            suffix="log.txt",
        ),
    group:
        "preproc"
    resources:
        mem_mb=config["total_mem_mb"],  #TODO update this, along with threads.. 
        disk_mb=1000000,  #1TB
    threads: 16
    script:
        "../scripts/bioformats_to_zarr.py"


rule concat_tiles:
    """ read in zarrs created for each tile, and write out as a single zarr"""
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
            datatype="concat_tiles",
            subject="{subject}",
            sample="{sample}",
            acq="{acq}",
            suffix="benchmark.tsv",
        )
    log:
        bids(
            root="logs",
            datatype="concat_tiles",
            subject="{subject}",
            sample="{sample}",
            acq="{acq}",
            suffix="log.txt",
        ),
    group:
        "preproc"
    resources:
        mem_mb=config["total_mem_mb"],
        disk_mb=1000000,  #1TB
    threads: 32
    container:
        None
    script:
        "../scripts/concat_tiles.py"
