def cmd_get_dataset(wildcards, input, output):
    cmds = []
    import tarfile

    # supports tar, tar.gz, tgz, zip, or folder name
    dataset_path = Path(input.dataset_path)
    suffix = dataset_path.suffix
    if dataset_path.is_dir():
        # we have a directory:
        # return command to copy folder
        cmds.append(f"ln -sr {input} {output}")

    elif tarfile.is_tarfile(dataset_path):
        # we have a tar file
        # check if gzipped:
        cmds.append(f"mkdir -p {output}")
        if suffix == "gz" or suffix == "tgz":
            cmds.append(f"tar -xzf {input} -C {output}")
        else:
            cmds.append(f"tar -xf {input} -C {output}")

    else:
        print(f"unsupported input: {dataset_path}")

    return " && ".join(cmds)


rule get_dataset:
    input:
        dataset_path=get_dataset_path,
    params:
        cmd=cmd_get_dataset,
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
                    suffix="spim",
                )
            )
        ),
    group:
        "preproc"
    shell:
        "{params.cmd}"


rule raw_to_metadata:
    input:
        ome_dir=rules.get_dataset.output.ome_dir,
    params:
        in_tif_glob=lambda wildcards, input: os.path.join(
            input.ome_dir,
            config["import"]["raw_tif_pattern"],
        ),
    output:
        metadata_json=bids(
            root=root,
            subject="{subject}",
            datatype="micr",
            sample="{sample}",
            acq="{acq}",
            suffix="spim.metadata.json",
        ),
    benchmark:
        bids(
            root="benchmarks",
            datatype="raw_to_metdata",
            subject="{subject}",
            sample="{sample}",
            acq="{acq}",
            suffix="benchmark.tsv",
        )
    group:
        "preproc"
    container:
        config["containers"]["spimprep"]
    script:
        "../scripts/raw_to_metadata.py"


rule tif_to_zarr:
    """ use dask to load tifs in parallel and write to zarr 
        output shape is (tiles,channels,z,y,x), with the 2d 
        images as the chunks"""
    input:
        ome_dir=rules.get_dataset.output.ome_dir,
        metadata_json=rules.raw_to_metadata.output.metadata_json,
    params:
        in_tif_glob=lambda wildcards, input: os.path.join(
            input.ome_dir,
            config["import"]["raw_tif_glob"],
        ),
        intensity_rescaling=config["import"]["intensity_rescaling"],
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
                    suffix="spim.zarr",
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
    group:
        "preproc"
    threads: 32
    container:
        config["containers"]["spimprep"]
    script:
        "../scripts/tif_to_zarr.py"
