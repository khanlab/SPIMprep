rule multiview_stitcher:
    input:
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
        registration_binning=config['multiview_stitcher']['registration_binning'],
        reg_channel_index=1, #later can make this a parameter chosen based on stains
    output:
        tiling_qc_png=bids(
                    root=root,
                    subject="{subject}",
                    datatype="qc",
                    sample="{sample}",
                    acq="{acq}",
                    desc="{desc}mvstitcher",
                    suffix="affinemetadata.tiles.png",
            ),
        ome_zarr=directory(
                bids(
                    root=root,
                    subject="{subject}",
                    datatype="micr",
                    sample="{sample}",
                    acq="{acq}",
                    desc="{desc}mvstitcher",
                    suffix="SPIM.ome.zarr",
                )
            )
    shadow: 'minimal'  #don't make this shadow, as we then have to copy files after (less efficient)
    #instead just make a work temp() folder for the intermediate zarrs - better for cloud too..
    benchmark:
        bids(
            root="benchmarks",
            datatype="multiview_stitcher",
            subject="{subject}",
            sample="{sample}",
            acq="{acq}",
            desc="{desc}",
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
            suffix="log.txt",
        ),
    threads: config["total_cores"]
    resources:
        mem_mb=config["total_mem_mb"],
    group:
        "preproc"
    script:
        "../scripts/multiview_stitcher.py"


