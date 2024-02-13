rule zarr_to_bdv:
    input:
        zarr=bids(
            root=work,
            subject="{subject}",
            datatype="micr",
            sample="{sample}",
            acq="{acq}",
            desc="{desc}",
            suffix="spim.zarr",
        ),
        metadata_json=rules.raw_to_metadata.output.metadata_json,
    params:
        max_downsampling_layers=5,  #1,2,4,8,16
        temp_h5=str(
            Path(
                bids(
                    root=work,
                    subject="{subject}",
                    datatype="micr",
                    sample="{sample}",
                    acq="{acq}",
                    desc="{desc}",
                    suffix="bdv",
                )
            )
            / "dataset.h5"
        ),
        #only temporary, is promptly deleted 
        temp_xml=str(
            Path(
                bids(
                    root=work,
                    subject="{subject}",
                    datatype="micr",
                    sample="{sample}",
                    acq="{acq}",
                    desc="{desc}",
                    suffix="bdv",
                )
            )
            / "dataset.xml"
        ),
        #only temporary, is promptly deleted 
    output:
        bdv_n5=temp(
            directory(
                bids(
                    root=work,
                    subject="{subject}",
                    datatype="micr",
                    sample="{sample}",
                    acq="{acq}",
                    desc="{desc}",
                    suffix="bdv.n5",
                )
            )
        ),
        bdv_xml=temp(
            bids(
                root=work,
                subject="{subject}",
                datatype="micr",
                sample="{sample}",
                acq="{acq}",
                desc="{desc}",
                suffix="bdv.xml",
            )
        ),
    benchmark:
        bids(
            root="benchmarks",
            datatype="zarr_to_n5_bdv",
            subject="{subject}",
            sample="{sample}",
            acq="{acq}",
            desc="{desc}",
            suffix="benchmark.tsv",
        )
    threads: 32
    group:
        "preproc"
    container:
        config["containers"]["spimprep"]
    script:
        "../scripts/zarr_to_n5_bdv.py"


def get_fiji_launcher_cmd(wildcards, output, threads, resources):
    launcher_opts_find = "-Xincgc"
    launcher_opts_replace = f"-XX:+UseG1GC -verbose:gc -XX:+PrintGCDateStamps -XX:ActiveProcessorCount={threads}"
    pipe_cmds = []
    pipe_cmds.append("ImageJ-linux64 --dry-run --headless --console")
    pipe_cmds.append(f"sed 's/{launcher_opts_find}/{launcher_opts_replace}'/")
    pipe_cmds.append(
        f"sed 's/-Xmx[0-9a-z]\+/-Xmx{resources.mem_mb}m -Xms{resources.mem_mb}m/'"
    )
    pipe_cmds.append("tr --delete '\\n'")
    return "|".join(pipe_cmds) + f" > {output.launcher} && chmod a+x {output.launcher} "


def get_macro_args_bigstitcher(wildcards, input, output):
    return "{dataset_xml} {ds_x} {ds_y} {ds_z} {min_r}".format(
        dataset_xml=input.dataset_xml,
        ds_x=config["bigstitcher"]["calc_pairwise_shifts"]["downsample_in_x"],
        ds_y=config["bigstitcher"]["calc_pairwise_shifts"]["downsample_in_y"],
        ds_z=config["bigstitcher"]["calc_pairwise_shifts"]["downsample_in_z"],
        min_r=config["bigstitcher"]["filter_pairwise_shifts"]["min_r"],
    )


rule bigstitcher:
    input:
        dataset_n5=rules.zarr_to_bdv.output.bdv_n5,
        dataset_xml=rules.zarr_to_bdv.output.bdv_xml,
        ijm=Path(workflow.basedir) / "macros" / "AutostitchMacro.ijm",
    params:
        fiji_launcher_cmd=get_fiji_launcher_cmd,
        macro_args=get_macro_args_bigstitcher,
        rm_old_xml=lambda wildcards, output: f"rm -f {output.dataset_xml}~?",
    output:
        launcher=temp(
            bids(
                root=work,
                subject="{subject}",
                datatype="micr",
                sample="{sample}",
                acq="{acq}",
                desc="{desc}",
                suffix="bigstitcherproc.sh",
            )
        ),
        dataset_xml=temp(
            bids(
                root=work,
                subject="{subject}",
                datatype="micr",
                sample="{sample}",
                acq="{acq}",
                desc="{desc}",
                suffix="bigstitcher.xml",
            )
        ),
    benchmark:
        bids(
            root="benchmarks",
            datatype="bigstitcherproc",
            subject="{subject}",
            sample="{sample}",
            acq="{acq}",
            desc="{desc}",
            suffix="benchmark.tsv",
        )
    log:
        bids(
            root="logs",
            datatype="bigstitcherproc",
            subject="{subject}",
            sample="{sample}",
            acq="{acq}",
            desc="{desc}",
            suffix="log.txt",
        ),
    container:
        config["containers"]["spimprep"]
    resources:
        runtime=30,  #this should be proportional to the number of images and image size
        mem_mb=10000,
    threads: 32
    group:
        "preproc"
    shell:
        "cp {input.dataset_xml} {output.dataset_xml} && "
        " {params.fiji_launcher_cmd} && "
        " echo ' -macro {input.ijm} \"{params.macro_args}\"' >> {output.launcher} "
        " && {output.launcher} &> {log} && {params.rm_old_xml}"


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


rule fuse_dataset:
    input:
        dataset_n5=bids(
            root=work,
            subject="{subject}",
            datatype="micr",
            sample="{sample}",
            acq="{acq}",
            desc="{desc}",
            suffix="bdv.n5",
        ),
        dataset_xml=bids(
            root=work,
            subject="{subject}",
            datatype="micr",
            sample="{sample}",
            acq="{acq}",
            desc="{desc}",
            suffix="bigstitcher.xml",
        ),
        ijm=Path(workflow.basedir) / "macros" / "FuseImageMacroZarr.ijm",
    params:
        fiji_launcher_cmd=get_fiji_launcher_cmd,
        macro_args=get_macro_args_zarr_fusion,
    output:
        launcher=temp(
            bids(
                root=work,
                subject="{subject}",
                datatype="micr",
                sample="{sample}",
                acq="{acq}",
                desc="stitched{desc}",
                stain="{stain}",
                suffix="fuseimagen5.sh",
            )
        ),
        zarr=temp(
            directory(
                bids(
                    root=work,
                    subject="{subject}",
                    datatype="micr",
                    sample="{sample}",
                    acq="{acq}",
                    desc="stitched{desc}",
                    stain="{stain}",
                    suffix="spim.zarr",
                )
            )
        ),
    benchmark:
        bids(
            root="benchmarks",
            datatype="fuse_dataset",
            subject="{subject}",
            sample="{sample}",
            acq="{acq}",
            desc="stitched{desc}",
            stain="{stain}",
            suffix="benchmark.tsv",
        )
    log:
        bids(
            root="logs",
            datatype="fuse_dataset",
            subject="{subject}",
            sample="{sample}",
            acq="{acq}",
            desc="stitched{desc}",
            stain="{stain}",
            suffix="log.txt",
        ),
    container:
        config["containers"]["spimprep"]
    resources:
        runtime=30,
        mem_mb=20000,
    threads: 32
    group:
        "preproc"
    shell:
        " {params.fiji_launcher_cmd} && "
        " echo ' -macro {input.ijm} \"{params.macro_args}\"' >> {output.launcher} "
        " && {output.launcher} &> {log}"


rule fuse_dataset_spark:
    input:
        dataset_n5=bids(
            root=work,
            subject="{subject}",
            datatype="micr",
            sample="{sample}",
            acq="{acq}",
            desc="{desc}",
            suffix="bdv.n5",
        ),
        dataset_xml=bids(
            root=work,
            subject="{subject}",
            datatype="micr",
            sample="{sample}",
            acq="{acq}",
            desc="{desc}",
            suffix="bigstitcher.xml",
        ),
        ijm=Path(workflow.basedir) / "macros" / "FuseImageMacroZarr.ijm",
    params:
    output:
        zarr=temp(
            directory(
                bids(
                    root=work,
                    subject="{subject}",
                    datatype="micr",
                    sample="{sample}",
                    acq="{acq}",
                    desc="sparkstitched{desc}",
                    stain="{stain}",
                    suffix="spim.zarr",
                )
            )
        ),
    benchmark:
        bids(
            root="benchmarks",
            datatype="fuse_dataset_spark",
            subject="{subject}",
            sample="{sample}",
            acq="{acq}",
            desc="stitched{desc}",
            stain="{stain}",
            suffix="benchmark.tsv",
        )
    log:
        bids(
            root="logs",
            datatype="fuse_dataset_spark",
            subject="{subject}",
            sample="{sample}",
            acq="{acq}",
            desc="stitched{desc}",
            stain="{stain}",
            suffix="log.txt",
        ),
    container:
        config["containers"]["spimprep"]
    resources:
        runtime=30,
        mem_mb=20000,
    threads: 32
    group:
        "preproc"
    shell:
        "affine-fusion ..."


