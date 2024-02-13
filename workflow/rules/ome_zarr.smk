rule zarr_to_ome_zarr:
    input:
        zarr=bids(
            root=work,
            subject="{subject}",
            datatype="micr",
            sample="{sample}",
            acq="{acq}",
            desc="{desc}",
            stain="{stain}",
            suffix="spim.zarr",
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
                    desc="{desc}",
                    stain="{stain}",
                    suffix="spim.ome.zarr",
                )
            )
        ),
    threads: 32
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
            root=work,
            subject="{subject}",
            datatype="micr",
            sample="{sample}",
            acq="{acq}",
            desc="{desc}",
            stain="{stain}",
            suffix="spim.ome.zarr.zip",
        ),
    output:
        nii=bids(
            root=root,
            subject="{subject}",
            datatype="micr",
            sample="{sample}",
            acq="{acq}",
            desc="{desc}",
            stain="{stain}",
            level="{level}",
            suffix="spim.nii",
        ),
    benchmark:
        bids(
            root="benchmarks",
            datatype="ome_zarr_to_nifti",
            subject="{subject}",
            sample="{sample}",
            acq="{acq}",
            desc="{desc}",
            stain="{stain}",
            level="{level}",
            suffix="benchmark.tsv",
        )
    group:
        "preproc"
    threads: 32
    container:
        config["containers"]["spimprep"]
    script:
        "../scripts/ome_zarr_to_nii.py"


rule zarr_masking_wip:
    input:
        zarr=rules.zarr_to_ome_zarr.output.zarr,
    output:
        zarr=directory(
            bids(
                root=root,
                subject="{subject}",
                datatype="micr",
                sample="{sample}",
                acq="{acq}",
                desc="{desc}",
                stain="{stain}",
                suffix="mask.ome.zarr",
            )
        ),
    container:
        config["containers"]["spimprep"]
    notebook:
        "../notebooks/zarr_masking.py.ipynb"
