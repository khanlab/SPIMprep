import tarfile


# targets
def get_all_targets():
    targets = []
    for i in range(len(datasets)):
        targets.extend(
            expand(
                bids(
                    root=root,
                    subject="{subject}",
                    datatype="micr",
                    sample="{sample}",
                    acq="{acq}",
                    suffix="SPIM.{extension}",
                ),
                subject=datasets.loc[i, "subject"],
                sample=datasets.loc[i, "sample"],
                acq=datasets.loc[i, "acq"],
                extension=(
                    "ome.zarr.zip" if config["ome_zarr"]["use_zipstore"] else "ome.zarr"
                ),
            )
        )
        targets.extend(
            expand(
                bids(
                    root=resampled,
                    subject="{subject}",
                    datatype="micr",
                    sample="{sample}",
                    acq="{acq}",
                    res="{level}x",
                    stain="{stain}",
                    suffix="SPIM.nii",
                ),
                subject=datasets.loc[i, "subject"],
                sample=datasets.loc[i, "sample"],
                acq=datasets.loc[i, "acq"],
                level=config["nifti"]["levels"],
                stain=[datasets.loc[i, "stain_0"], datasets.loc[i, "stain_1"]],
            )
        )

    return targets


def get_bids_toplevel_targets():
    targets = []
    targets.append(Path(root) / "README.md")
    targets.append(Path(root) / "dataset_description.json")
    targets.append(Path(root) / "samples.tsv")
    targets.append(Path(root) / "samples.json")
    targets.append(Path(resampled) / "dataset_description.json")
    return targets


def get_input_dataset(wildcards):
    """returns path to extracted dataset or path to provided input folder"""

    in_dataset_path = get_dataset_path(wildcards)

    ## this logic below became messy when trying to deal with both local and remote files
    # can simplify perhaps if I use cloudpathlib (or some way to check if is dir or tar remotely)

    # first check if we have a google cloud storage URI, and split path
    split_path = in_dataset_path.split("://")
    if len(split_path) > 1:
        # we have a storage uri, so we need to check if it is a tarfile
        # or a directory
        dataset_path = Path(split_path[1])

        suffix = dataset_path.suffix

        if dataset_path.suffix == ".tar" or dataset_path.suffix == ".tgz":
            # we have a directory already, just point to it
            print("is a tar")
            return rules.extract_dataset.output.ome_dir.format(**wildcards)

        elif tarfile.is_tarfile(dataset_path):
            # dataset was a tar file, so point to the extracted folder
            print("is a dir")
            return str(dataset_path)

        else:
            print(f"unsupported input: {dataset_path}")

    else:
        dataset_path = Path(in_dataset_path)

        suffix = dataset_path.suffix

        if dataset_path.is_dir():
            # we have a directory already, just point to it
            print("is a dir")
            return str(dataset_path)

        elif tarfile.is_tarfile(dataset_path):
            print("is a tar")
            # dataset was a tar file, so point to the extracted folder
            return rules.extract_dataset.output.ome_dir.format(**wildcards)

        else:
            print(f"unsupported input: {dataset_path}")


# import
def cmd_extract_dataset(wildcards, input, output):
    cmds = []

    # supports tar, tar.gz, tgz, or folder name
    dataset_path = Path(input.dataset_path)
    suffix = dataset_path.suffix
    if dataset_path.is_dir():
        # we have a directory
        print("input directory not copied/extracted by this rule")

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


def get_dataset_path(wildcards):
    df = datasets.query(
        f"subject=='{wildcards.subject}' and sample=='{wildcards.sample}' and acq=='{wildcards.acq}'"
    )
    return df.dataset_path.to_list()[0]


def get_dataset_storage_path(wildcards):
    dataset_path = get_dataset_path(wildcards)
    split_path = dataset_path.split("://")
    if len(split_path) > 1:
        return storage(dataset_path)
    else:
        return dataset_path


def get_stains(wildcards):
    df = datasets.query(
        f"subject=='{wildcards.subject}' and sample=='{wildcards.sample}' and acq=='{wildcards.acq}'"
    )

    return [
        df.stain_0.to_list()[0],
        df.stain_1.to_list()[0],
    ]


# bids
def bids_tpl(root, template, **entities):
    """bids() wrapper for files in tpl-template folder"""
    return str(Path(bids(root=root, tpl=template)) / bids(tpl=template, **entities))


# bigstitcher
def get_fiji_launcher_cmd(wildcards, output, threads, resources):
    launcher_opts_find = "-Xincgc"
    launcher_opts_replace = f"-XX:+UseG1GC -verbose:gc -XX:+PrintGCDateStamps -XX:ActiveProcessorCount={threads}"
    pattern_mem = r"-Xmx[0-9a-z]\+"
    pipe_cmds = []
    pipe_cmds.append("ImageJ-linux64 --dry-run --headless --console")
    pipe_cmds.append(f"sed 's/{launcher_opts_find}/{launcher_opts_replace}'/")
    pipe_cmds.append(
        f"sed 's/{pattern_mem}/-Xmx{resources.mem_mb}m -Xms{resources.mem_mb}m/'"
    )
    pipe_cmds.append("tr --delete '\\n'")
    return "|".join(pipe_cmds) + f" > {output.launcher} && chmod a+x {output.launcher} "


def get_macro_args_bigstitcher(wildcards, input, output):
    return "{dataset_xml} {ds_x} {ds_y} {ds_z} {min_r}".format(
        dataset_xml=output.dataset_xml,
        ds_x=config["bigstitcher"]["calc_pairwise_shifts"]["downsample_in_x"],
        ds_y=config["bigstitcher"]["calc_pairwise_shifts"]["downsample_in_y"],
        ds_z=config["bigstitcher"]["calc_pairwise_shifts"]["downsample_in_z"],
        min_r=config["bigstitcher"]["filter_pairwise_shifts"]["min_r"],
    )


def get_macro_args_zarr_fusion(wildcards, input, output):
    return "{dataset_xml} {downsampling} {channel:02d} {output_zarr} {bsx} {bsy} {bsz} {bsfx} {bsfy} {bsfz}".format(
        dataset_xml=input.dataset_xml,
        downsampling=config["bigstitcher"]["fuse_dataset"]["downsampling"],
        channel=get_stains(wildcards).index(wildcards.stain),
        output_zarr=output.zarr,
        bsx=config["bigstitcher"]["fuse_dataset"]["block_size_x"],
        bsy=config["bigstitcher"]["fuse_dataset"]["block_size_y"],
        bsz=config["bigstitcher"]["fuse_dataset"]["block_size_z"],
        bsfx=config["bigstitcher"]["fuse_dataset"]["block_size_factor_x"],
        bsfy=config["bigstitcher"]["fuse_dataset"]["block_size_factor_y"],
        bsfz=config["bigstitcher"]["fuse_dataset"]["block_size_factor_z"],
    )
