import tarfile
from snakebids import bids as _bids
from upath import UPath as Path
from lib.cloud_io import is_remote, is_remote_gcs


def bids(root, *args, **kwargs):
    """replacement to _bids() that adds the storage() tag
    if the path points is remote (eg s3, gcs)"""
    bids_str = _bids(root=root, *args, **kwargs)

    if is_remote(root):
        return storage(bids_str)
    else:
        return bids_str


def expand_bids(expand_kwargs, **bids_kwargs):
    """replacement to expand(_bids()) that adds the storage() tag
    if the path points is remote (eg s3, gcs), **after**
    expanding (since it is not allowed before, hence why this
    function is needed)"""

    files = expand(_bids(**bids_kwargs), **expand_kwargs)

    if is_remote(files[0]):
        return [storage(f) for f in files]
    else:
        return files


def remote_file(filename):
    if is_remote(filename):
        return storage(str(filename))
    else:
        return filename


def remote_directory(dirname):
    if is_remote(dirname):
        return storage(directory(str(dirname)))
    else:
        return directory(dirname)


def directory_bids(root, *args, **kwargs):
    """Similar to expand_bids, this replacement function
    is needed to ensure storage() comes after directory() tags"""
    bids_str = _bids(root=root, *args, **kwargs)

    if is_remote(root):
        return storage(directory(bids_str))
    else:
        return directory(bids_str)


def bids_toplevel(root, filename):
    """This obtains a path for a file at the top-level of the bids
    dataset, applying the storage() tag if it is remote"""

    bids_str = str(Path(_bids(root=root)) / filename)

    if is_remote(root):
        return storage(bids_str)
    else:
        return bids_str


def get_extension_ome_zarr():
    """This function returns a ome.zarr extension. If the file is
    remote, it appends a .snakemake_touch to the path so that the rule
     can process the file remotely without having snakemake copy it over
    (e.g. the touch file is used instead). Also appends a .zip if the
    zipstore option is enabled."""

    if is_remote(config["root"]):
        return "ome.zarr/.snakemake_touch"
    else:
        if config["ome_zarr"]["use_zipstore"]:
            return "ome.zarr.zip"
        else:
            return "ome.zarr"


# targets
def get_all_targets():
    targets = []
    for i in range(len(datasets)):
        targets.extend(
            expand_bids(
                root=root,
                subject="{subject}",
                datatype="micr",
                sample="{sample}",
                acq="{acq}",
                suffix="SPIM.{extension}",
                expand_kwargs=dict(
                    subject=datasets.loc[i, "subject"],
                    sample=datasets.loc[i, "sample"],
                    acq=datasets.loc[i, "acq"],
                    extension=[get_extension_ome_zarr(), "json"],
                ),
            )
        )
        targets.extend(
            expand_bids(
                root=resampled,
                subject="{subject}",
                datatype="micr",
                sample="{sample}",
                acq="{acq}",
                res="{level}x",
                stain="{stain}",
                suffix="SPIM.nii",
                expand_kwargs=dict(
                    subject=datasets.loc[i, "subject"],
                    sample=datasets.loc[i, "sample"],
                    acq=datasets.loc[i, "acq"],
                    level=config["nifti"]["levels"],
                    stain=get_stains_by_row(i),
                ),
            )
        )
    return targets


def get_all_subj_html(wildcards):
    htmls = []

    for i in range(len(datasets)):

        html = "{root}/qc/sub-{subject}_sample-{sample}_acq-{acq}/subject.html".format(
            root=root,
            subject=datasets.loc[i, "subject"],
            sample=datasets.loc[i, "sample"],
            acq=datasets.loc[i, "acq"],
        )
        htmls.append(remote_file(html))

    return htmls


def get_bids_toplevel_targets():
    targets = []
    targets.append(bids_toplevel(root, "README.md"))
    targets.append(bids_toplevel(root, "dataset_description.json"))
    targets.append(bids_toplevel(root, "samples.tsv"))
    targets.append(bids_toplevel(root, "samples.json"))
    targets.append(bids_toplevel(resampled, "dataset_description.json"))
    return targets


def get_qc_targets():
    targets = []
    if config["report"]["create_report"]:
        targets.append(remote_file(Path(root) / "qc" / "qc_report.html"))
        targets.append(remote_file(Path(root) / "qc" / "README.md"))
    return targets


def dataset_is_remote(wildcards):
    return is_remote_gcs(Path(get_dataset_path(wildcards)))


def get_input_dataset(wildcards):
    """returns path to extracted dataset or path to provided input folder"""
    dataset_path = Path(get_dataset_path(wildcards))
    suffix = dataset_path.suffix

    if is_remote_gcs(dataset_path):
        return rules.cp_from_gcs.output.ome_dir.format(**wildcards)

    if dataset_path.is_dir():
        return get_dataset_path_remote(wildcards)

    elif tarfile.is_tarfile(dataset_path):
        # dataset was a tar file, so point to the extracted folder
        return rules.extract_dataset.output.ome_dir.format(**wildcards)

    else:
        print(f"unsupported input: {dataset_path}")


def get_metadata_json(wildcards):
    """returns path to metadata, extracted from local or gcs"""
    dataset_path = Path(get_dataset_path(wildcards))
    suffix = dataset_path.suffix

    if is_remote_gcs(dataset_path):
        return rules.blaze_to_metadata_gcs.output.metadata_json.format(**wildcards)
    else:
        return rules.blaze_to_metadata.output.metadata_json.format(**wildcards)


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


def get_dataset_path_remote(wildcards):
    path = get_dataset_path(wildcards)
    if is_remote(path):
        return storage(path)
    else:
        return path


def get_dataset_path_gs(wildcards):
    path = Path(get_dataset_path(wildcards)).path
    return f"gs://{path}"


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
    launcher_opts_replace = f"-XX:+UseG1GC -verbose:gc -XX:+PrintGCDateStamps -XX:ActiveProcessorCount={threads} -XX:ParallelGCThreads={int(threads/2)} -XX:ConcGCThreads={int(threads/4)}"
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


def get_output_ome_zarr_uri():
    if is_remote(config["root"]):
        return _bids(
            root=root,
            subject="{subject}",
            datatype="micr",
            sample="{sample}",
            acq="{acq}",
            suffix="SPIM.ome.zarr",
        )
    else:
        return "local://" + _bids(
            root=root,
            subject="{subject}",
            datatype="micr",
            sample="{sample}",
            acq="{acq}",
            suffix="SPIM.ome.zarr",
        )


def get_output_ome_zarr(acq_type):
    if is_remote(config["root"]):
        return {
            "zarr": touch(
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
        }
    else:
        if config["write_ome_zarr_direct"]:
            return {
                "zarr": directory_bids(
                    root=root,
                    subject="{subject}",
                    datatype="micr",
                    sample="{sample}",
                    acq=f"{{acq,[a-zA-Z0-9]*{acq_type}[a-zA-Z0-9]*}}",
                    suffix="SPIM.ome.zarr",
                )
            }
        else:
            return {
                "zarr": temp(
                    directory_bids(
                        root=work,
                        subject="{subject}",
                        datatype="micr",
                        sample="{sample}",
                        acq=f"{{acq,[a-zA-Z0-9]*{acq_type}[a-zA-Z0-9]*}}",
                        suffix="SPIM.ome.zarr",
                    )
                )
            }


def get_input_ome_zarr_to_nii():
    if is_remote(config["root"]):
        return bids(
            root=root,
            subject="{subject}",
            datatype="micr",
            sample="{sample}",
            acq="{acq}",
            suffix="SPIM.{extension}".format(extension=get_extension_ome_zarr()),
        )
    else:
        if config["write_ome_zarr_direct"]:
            return bids(
                root=root,
                subject="{subject}",
                datatype="micr",
                sample="{sample}",
                acq="{acq}",
                suffix="SPIM.{extension}".format(extension=get_extension_ome_zarr()),
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
    protocol = Path(config["root"]).protocol
    if protocol == "gcs":
        # currently only works with gcs
        creds = os.path.expanduser(config["remote_creds"])
        return {"creds": creds}
    else:
        return {}
