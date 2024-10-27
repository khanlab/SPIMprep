rule zarr_to_bdv:
    input:
        zarr=bids(
            root=work,
            subject="{subject}",
            datatype="micr",
            sample="{sample}",
            acq="{acq}",
            desc="{desc}",
            suffix="SPIM.{ext}".format(
                ext="zarr.zip" if config["use_zipstore"] else "zarr"
            ),
        ),
        metadata_json=rules.copy_blaze_metadata.output.metadata_json,
    params:
        max_downsampling_layers=5,
        temp_h5=str(
            Path(
                bids(
                    root=workn5,
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
                    root=workn5,
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
                    root=workn5,
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
                root=workn5,
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
    threads: config["total_cores"]
    group:
        "preproc"
    container:
        config["containers"]["spimprep"]
    script:
        "../scripts/zarr_to_n5_bdv.py"


rule bigstitcher_stitching:
    input:
        dataset_n5=rules.zarr_to_bdv.output.bdv_n5,
        dataset_xml=rules.zarr_to_bdv.output.bdv_xml,
    params:
        downsampling="--downsampling={dsx},{dsy},{dsz}".format(
            dsx=config["bigstitcher"]["calc_pairwise_shifts"]["downsample_in_x"],
            dsy=config["bigstitcher"]["calc_pairwise_shifts"]["downsample_in_y"],
            dsz=config["bigstitcher"]["calc_pairwise_shifts"]["downsample_in_z"],
        ),
        min_r="--minR={min_r}".format(
            min_r=config["bigstitcher"]["filter_pairwise_shifts"]["min_r"]
        ),
        max_shift="--maxShiftTotal={max_shift}".format(
            max_shift=config["bigstitcher"]["filter_pairwise_shifts"][
                "max_shift_total"
            ]
        ),
        mem_gb=lambda wildcards, resources: "{mem_gb}".format(
            mem_gb=int(resources.mem_mb / 1000)
        ),
        rm_old_xml=lambda wildcards, output: f"rm -f {output.dataset_xml}~?",
    output:
        dataset_xml=temp(
            bids(
                root=workn5,
                subject="{subject}",
                datatype="micr",
                sample="{sample}",
                acq="{acq}",
                desc="{desc}",
                suffix="bigstitcherstitching.xml",
            )
        ),
    benchmark:
        bids(
            root="benchmarks",
            datatype="bigstitcherstitching",
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
        mem_mb=int(config["total_mem_mb"] * 0.9),
    threads: config["total_cores"]
    group:
        "preproc"
    shell:
        "cp {input.dataset_xml} {output.dataset_xml} && "
        "stitching {params.mem_gb} {threads} -x {output.dataset_xml} "
        " {params.min_r} {params.downsampling} {params.max_shift} && "
        "{params.rm_old_xml}"


rule bigstitcher_detect_interestpoints:
    input:
        dataset_n5=rules.zarr_to_bdv.output.bdv_n5,
        dataset_xml=rules.zarr_to_bdv.output.bdv_xml,
    params:
        downsampling="--downsampleXY={dsxy}".format(
            dsxy=config["bigstitcher"]["interest_points"]["downsample_xy"],
        ),
        min_intensity="--minIntensity={min_intensity}".format(
            min_intensity=config["bigstitcher"]["interest_points"]["min_intensity"]
        ),
        max_intensity="--maxIntensity={max_intensity}".format(
            max_intensity=config["bigstitcher"]["interest_points"]["max_intensity"]
        ),
        label="--label={label}".format(
            label=config["bigstitcher"]["interest_points"]["label"]
        ),
        threshold="--threshold={threshold}".format(
            threshold=config["bigstitcher"]["interest_points"]["threshold"]
        ),
        sigma="--sigma={sigma}".format(
            sigma=config["bigstitcher"]["interest_points"]["sigma"]
        ),
        mem_gb=lambda wildcards, resources: "{mem_gb}".format(
            mem_gb=int(resources.mem_mb / 1000)
        ),
        rm_old_xml=lambda wildcards, output: f"rm -f {output.dataset_xml}~?",
    output:
        dataset_xml=temp(
            bids(
                root=workn5,
                subject="{subject}",
                datatype="micr",
                sample="{sample}",
                acq="{acq}",
                desc="{desc}",
                suffix="bigstitcherdetectinterestpoints.xml",
            )
        ),
    benchmark:
        bids(
            root="benchmarks",
            datatype="bigstitcherdetectinterestpoints",
            subject="{subject}",
            sample="{sample}",
            acq="{acq}",
            desc="{desc}",
            suffix="benchmark.tsv",
        )
    log:
        bids(
            root="logs",
            datatype="bigstitcherdetectinterestpoints",
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
        mem_mb=int(config["total_mem_mb"] * 0.9),
    threads: config["total_cores"]
    group:
        "preproc"
    shell:
        "cp {input.dataset_xml} {output.dataset_xml} && "
        "detect-interestpoints {params.mem_gb} {threads} -x {output.dataset_xml} "
        " {params.min_intensity} {params.max_intensity} {params.downsampling} "
        " {params.sigma} {params.threshold} {params.label}  && "
        "{params.rm_old_xml}"


rule bigstitcher_match_interestpoints:
    input:
        dataset_n5=rules.zarr_to_bdv.output.bdv_n5,
        dataset_xml=rules.bigstitcher_detect_interestpoints.output.dataset_xml,
    params:
        label="--label={label}".format(label="beads"),
        method="--method={method}".format(method="ICP"),
        mem_gb=lambda wildcards, resources: "{mem_gb}".format(
            mem_gb=int(resources.mem_mb / 1000)
        ),
        rm_old_xml=lambda wildcards, output: f"rm -f {output.dataset_xml}~?",
    output:
        dataset_xml=temp(
            bids(
                root=workn5,
                subject="{subject}",
                datatype="micr",
                sample="{sample}",
                acq="{acq}",
                desc="{desc}",
                suffix="bigstitchermatchinterestpoints.xml",
            )
        ),
    benchmark:
        bids(
            root="benchmarks",
            datatype="bigstitchermatchinterestpoints",
            subject="{subject}",
            sample="{sample}",
            acq="{acq}",
            desc="{desc}",
            suffix="benchmark.tsv",
        )
    log:
        bids(
            root="logs",
            datatype="bigstitchermatchinterestpoints",
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
        mem_mb=int(config["total_mem_mb"] * 0.9),
    threads: config["total_cores"]
    group:
        "preproc"
    shell:
        "cp {input.dataset_xml} {output.dataset_xml} && "
        "match-interestpoints {params.mem_gb} {threads} -x {output.dataset_xml} "
        " {params.method} {params.label}  && "
        "{params.rm_old_xml}"


rule bigstitcher_solver:
    input:
        dataset_n5=rules.zarr_to_bdv.output.bdv_n5,
        dataset_xml=(
            rules.bigstitcher_match_interestpoints.output.dataset_xml
            if config["bigstitcher"]["use_interestpoints"]
            else rules.bigstitcher_stitching.output.dataset_xml
        ),
    params:
        label="--label={label}".format(label="beads"),
        lambdareg="--lambda={lambdareg}".format(lambdareg=0.1),
        source="-s {source}".format(
            source="IP" if config["bigstitcher"]["use_interestpoints"] else "STITCHING"
        ),
        method="--method={method}".format(
            method=config["bigstitcher"]["global_optimization"]["method"]
        ),
        mem_gb=lambda wildcards, resources: "{mem_gb}".format(
            mem_gb=int(resources.mem_mb / 1000)
        ),
        rm_old_xml=lambda wildcards, output: f"rm -f {output.dataset_xml}~?",
    output:
        dataset_xml=temp(
            bids(
                root=workn5,
                subject="{subject}",
                datatype="micr",
                sample="{sample}",
                acq="{acq}",
                desc="{desc}",
                suffix="bigstitchersolver.xml",
            )
        ),
    benchmark:
        bids(
            root="benchmarks",
            datatype="bigstitchersolver",
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
        mem_mb=int(config["total_mem_mb"] * 0.9),
    threads: config["total_cores"]
    group:
        "preproc"
    shell:
        "cp {input.dataset_xml} {output.dataset_xml} && "
        "solver {params.mem_gb} {threads} -x {output.dataset_xml} "
        " {params.source} {params.label} {params.lambdareg}"
        " {params.method} && "
        "{params.rm_old_xml}"


rule bigstitcher_fusion:
    input:
        dataset_n5=bids(
            root=workn5,
            subject="{subject}",
            datatype="micr",
            sample="{sample}",
            acq="{acq}",
            desc="{desc}",
            suffix="bdv.n5",
        ),
        dataset_xml=bids(
            root=workn5,
            subject="{subject}",
            datatype="micr",
            sample="{sample}",
            acq="{acq}",
            desc="{desc}",
            suffix="bigstitcher{}.xml".format(
                "solver"
                if config["bigstitcher"]["global_optimization"]["enabled"]
                else "stitching"
            ),
        ),
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
        mem_gb=lambda wildcards, resources: "{mem_gb}".format(
            mem_gb=int(resources.mem_mb / 1000)
        ),
    output:
        zarr=temp(
            directory(
                bids(
                    root=workn5,
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
        mem_mb=int(config["total_mem_mb"] * 0.9),
    threads: config["total_cores"]
    group:
        "preproc"
    shell:
        "affine-fusion {params.mem_gb} {threads} --preserveAnisotropy -x {input.dataset_xml} "
        " -o {output.zarr} -d /fused/s0 -s ZARR "
        " --UINT16 --minIntensity 0 --maxIntensity 65535 "
        "{params.block_size} {params.block_size_factor} {params.channel}"
