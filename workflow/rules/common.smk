import tarfile


def get_extension_ome_zarr():
    if config["write_to_remote"]:
        return "ome.zarr/.snakemake_touch"
    else:
        if config["ome_zarr"]["use_zipstore"]:
            return "ome.zarr.zip"
        else:
            return "ome.zarr"


def final(path_or_paths):
    if config["write_to_remote"]:
        if type(path_or_paths) == list:
            out_paths = []
            for path in path_or_paths:
                out_paths.append(storage(os.path.join(config["remote_prefix"], path)))
            return out_paths
        else:
            return storage(os.path.join(config["remote_prefix"], path_or_paths))
    else:
        return path_or_paths


# targets
def get_all_targets():
    targets = []
    for i in range(len(datasets)):
        targets.extend(
            final(
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
                    extension=[get_extension_ome_zarr(), "json"],
                )
            )
        )
        targets.extend(

            final(
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
                stain=get_stains_by_row(i),
            ))
        )

    return targets


def get_bids_toplevel_targets():
    targets = []
    targets.append(Path(root) / "README.md")
    targets.append(Path(root) / "dataset_description.json")
    targets.append(Path(root) / "samples.tsv")
    targets.append(Path(root) / "samples.json")
    targets.append(Path(resampled) / "dataset_description.json")
    return [final(target) for target in targets]


def get_input_dataset(wildcards):
    """returns path to extracted dataset or path to provided input folder"""
    in_dataset = get_dataset_path(wildcards)

    dataset_path = Path(get_dataset_path(wildcards))
    suffix = dataset_path.suffix

    if dataset_path.is_dir():
        # we have a directory already, just point to it
        return str(dataset_path)

    elif tarfile.is_tarfile(dataset_path):
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


def get_stains_by_row(i):

    # Select columns that match the pattern 'stain_'
    stain_columns = datasets.filter(like="stain_").columns

    # Select values for a given row
    return datasets.loc[i, stain_columns].dropna().tolist()


def get_stains(wildcards):
    df = datasets.query(
        f"subject=='{wildcards.subject}' and sample=='{wildcards.sample}' and acq=='{wildcards.acq}'"
    )

    # Select columns that match the pattern 'stain_'
    stain_columns = df.filter(like="stain_").columns

    return df.iloc[0][stain_columns].dropna().tolist()


# bigstitcher
def get_fiji_launcher_cmd(wildcards, output, threads, resources):
    launcher_opts_find = "-Xincgc"
    launcher_opts_replace = f"-XX:+UseG1GC -verbose:gc -XX:+PrintGCDateStamps -XX:ActiveProcessorCount={threads}"
    pipe_cmds = []
    pipe_cmds.append("ImageJ-linux64 --dry-run --headless --console")
    pipe_cmds.append(f"sed 's/{launcher_opts_find}/{launcher_opts_replace}'/")
    pipe_cmds.append(
        rf"sed 's/-Xmx[0-9a-z]\+/-Xmx{resources.mem_mb}m -Xms{resources.mem_mb}m/'"
    )
    pipe_cmds.append("tr --delete '\\n'")
    return "|".join(pipe_cmds) + f" > {output.launcher} && chmod a+x {output.launcher} "


def get_macro_args_bigstitcher(wildcards, input, output):
    return "{dataset_xml} {pairwise_method} {ds_x} {ds_y} {ds_z} {do_filter} {min_r} {do_global} {global_strategy}".format(
        dataset_xml=output.dataset_xml,
        pairwise_method=config["bigstitcher"]["calc_pairwise_shifts"]["methods"][
            config["bigstitcher"]["calc_pairwise_shifts"]["method"]
        ],
        ds_x=config["bigstitcher"]["calc_pairwise_shifts"]["downsample_in_x"],
        ds_y=config["bigstitcher"]["calc_pairwise_shifts"]["downsample_in_y"],
        ds_z=config["bigstitcher"]["calc_pairwise_shifts"]["downsample_in_z"],
        do_filter=config["bigstitcher"]["filter_pairwise_shifts"]["enabled"],
        min_r=config["bigstitcher"]["filter_pairwise_shifts"]["min_r"],
        do_global=config["bigstitcher"]["global_optimization"]["enabled"],
        global_strategy=config["bigstitcher"]["global_optimization"]["strategies"][
            config["bigstitcher"]["global_optimization"]["strategy"]
        ],
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


def get_output_ome_zarr(acq_type):
    if config["write_to_remote"]:
        return {
            "zarr": touch(
                final(
                    bids(
                        root=root,
                        subject="{subject}",
                        datatype="micr",
                        sample="{sample}",
                        acq=f"{{acq,[a-zA-Z0-9]*{acq_type}[a-zA-Z0-9]*}}",
                        suffix="SPIM.{extension}".format(
                            extension=get_extension_ome_zarr()
                        ),
                    )
                )
            )
        }
    else:
        if config["write_ome_zarr_direct"]:
            return {
                "zarr": final(
                    directory(
                        bids(
                            root=root,
                            subject="{subject}",
                            datatype="micr",
                            sample="{sample}",
                            acq=f"{{acq,[a-zA-Z0-9]*{acq_type}[a-zA-Z0-9]*}}",
                            suffix="SPIM.ome.zarr",
                        )
                    )
                )
            }
        else:
            return {
                "zarr": temp(
                    directory(
                        bids(
                            root=work,
                            subject="{subject}",
                            datatype="micr",
                            sample="{sample}",
                            acq=f"{{acq,[a-zA-Z0-9]*{acq_type}[a-zA-Z0-9]*}}",
                            suffix="SPIM.ome.zarr",
                        )
                    )
                )
            }


def get_input_ome_zarr_to_nii():
    if config["write_to_remote"]:
        return final(
            bids(
                root=root,
                subject="{subject}",
                datatype="micr",
                sample="{sample}",
                acq="{acq}",
                suffix="SPIM.{extension}".format(extension=get_extension_ome_zarr()),
            )
        )
    else:
        if config["write_ome_zarr_direct"]:
            return final(
                bids(
                    root=root,
                    subject="{subject}",
                    datatype="micr",
                    sample="{sample}",
                    acq="{acq}",
                    suffix="SPIM.{extension}".format(
                        extension=get_extension_ome_zarr()
                    ),
                )
            )
        else:
            return bids(
                root=work,
                subject="{subject}",
                datatype="micr",
                sample="{sample}",
                acq="{acq}",
                suffix=f"SPIM.ome.zarr",
            )


def get_storage_creds():
    """for rules that deal with remote storage directly"""
    if config["write_to_remote"]:
        # currently only works with gcs
        creds = os.path.expanduser(config["remote_creds"])
        return {"creds": creds}
    else:
        return {}
