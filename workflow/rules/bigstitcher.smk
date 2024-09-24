rule zarr_to_bdv:
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
        max_downsampling_layers=5,
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
    log:
        bids(
            root="logs",
            datatype="zarr_to_n5_bdv",
            subject="{subject}",
            sample="{sample}",
            acq="{acq}",
            desc="{desc}",
            suffix="log.txt",
        ),
    threads: config["cores_per_rule"]
    group:
        "preproc"
    container:
        config["containers"]["spimprep"]
    script:
        "../scripts/zarr_to_n5_bdv.py"


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
        runtime=30,
        mem_mb=10000,
    threads: config["cores_per_rule"]
    group:
        "preproc"
    shell:
        "cp {input.dataset_xml} {output.dataset_xml} && "
        " {params.fiji_launcher_cmd} && "
        " echo ' -macro {input.ijm} \"{params.macro_args}\"' >> {output.launcher} "
        " && {output.launcher} |& tee {log} && {params.rm_old_xml}"


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
                    suffix="SPIM.zarr",
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
    threads: config["cores_per_rule"]
    group:
        "preproc"
    shell:
        " {params.fiji_launcher_cmd} && "
        " echo ' -macro {input.ijm} \"{params.macro_args}\"' >> {output.launcher} "
        " && {output.launcher} |& tee {log}"


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
        channel=lambda wildcards: "--channelId={channel}".format(
            channel=get_stains(wildcards).index(wildcards.stain)
        ),
        block_size="--blockSize={bsx},{bsy},{bsz}".format(
            bsx=config["bigstitcher"]["fuse_dataset"]["block_size_x"],
            bsy=config["bigstitcher"]["fuse_dataset"]["block_size_y"],
            bsz=config["bigstitcher"]["fuse_dataset"]["block_size_z"],
        ),
        block_size_factor="--blockScale={bsfx},{bsfy},{bsfz}".format(
            bsfx=config["bigstitcher"]["fuse_dataset"]["block_size_factor_x"],
            bsfy=config["bigstitcher"]["fuse_dataset"]["block_size_factor_y"],
            bsfz=config["bigstitcher"]["fuse_dataset"]["block_size_factor_z"],
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
                    desc="sparkstitched{desc}",
                    stain="{stain}",
                    suffix="SPIM.zarr",
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
    threads: config["cores_per_rule"]
    group:
        "preproc"
    shell:
        "affine-fusion --preserveAnisotropy -x {input.dataset_xml} "
        " -o {output.zarr} -d /fused/s0 -s ZARR "
        " --UINT16 --minIntensity 0 --maxIntensity 65535 "


#        " --UINT8 --minIntensity 0 --maxIntensity 255 "
"{params}"  # all the params
