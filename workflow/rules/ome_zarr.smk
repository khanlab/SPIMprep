rule zarr_to_ome_zarr:
    input:
        zarr=lambda wildcards: expand(
            bids(
                root=work,
                subject="{subject}",
                datatype="micr",
                sample="{sample}",
                acq="{acq}",
                desc="{desc}",
                stain="{stain}",
                suffix="spim.zarr",
            ),
            stain=get_stains(wildcards),
            desc=config["ome_zarr"]["desc"],
            allow_missing=True,
        ),
        metadata_json=rules.raw_to_metadata.output.metadata_json,
    params:
        max_downsampling_layers=config["ome_zarr"]["max_downsampling_layers"],
        rechunk_size=config["ome_zarr"]["rechunk_size"],
        scaling_method=config["ome_zarr"]["scaling_method"],
        downsampling=config["bigstitcher"]["fuse_dataset"]["downsampling"],
    output:
        zarr=temp(
            directory(
                bids(
                    root=work,
                    subject="{subject}",
                    datatype="micr",
                    sample="{sample}",
                    acq="{acq}",
                    suffix="spim.ome.zarr",
                )
            )
        ),
    threads: 32
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


rule ome_zarr_to_zipstore:
    """ generic rule to process any ome.zarr from work """
    input:
        zarr=f"{work}/{{prefix}}.ome.zarr",
    output:
        zarr_zip=f"{root}/{{prefix}}.ome.zarr.zip",
    log:
        "logs/ome_zarr_to_zipstore/{prefix}.log",
    group:
        "preproc"
    shell:
        "7z a -mx0 -tzip {output.zarr_zip} {input.zarr}/. &> {log}"


rule ome_zarr_to_nii:
    input:
        zarr=bids(
            root=root,
            subject="{subject}",
            datatype="micr",
            sample="{sample}",
            acq="{acq}",
            suffix="spim.ome.zarr.zip",
        ),
    params:
        channel_index=lambda wildcards: get_stains(wildcards).index(wildcards.stain),
    output:
        nii=bids(
            root=resampled,
            subject="{subject}",
            datatype="micr",
            sample="{sample}",
            acq="{acq}",
            res="{level}x",
            stain="{stain}",
            suffix="spim.nii",
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
    threads: 32
    container:
        config["containers"]["spimprep"]
    script:
        "../scripts/ome_zarr_to_nii.py"
