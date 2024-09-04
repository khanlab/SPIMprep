
rule fit_basic_flatfield_corr:
    """ BaSiC flatfield correction"""
    input:
        zarr=lambda wildcards: bids(
            root=work,
            subject="{subject}",
            datatype="micr",
            sample="{sample}",
            acq="{acq}",
            desc="rawfromgcs" if dataset_is_remote(wildcards) else "raw",
            suffix="SPIM.zarr",
        ).format(**wildcards),
    params:
        channel=lambda wildcards: get_stains(wildcards).index(wildcards.stain),
        max_n_images=config["basic_flatfield_corr"]["max_n_images"],
        basic_opts=config["basic_flatfield_corr"]["fitting_opts"],
    output:
        model_dir=temp(
            directory(
                bids(
                    root=work,
                    subject="{subject}",
                    datatype="micr",
                    sample="{sample}",
                    acq="{acq}",
                    stain="{stain}",
                    suffix="basicmodel",
                )
            )
        ),
    resources:
        runtime=90,
        mem_mb=64000,
    threads: 8
    benchmark:
        bids(
            root="benchmarks",
            datatype="fit_basic_flatfield",
            subject="{subject}",
            sample="{sample}",
            acq="{acq}",
            stain="{stain}",
            suffix="benchmark.tsv",
        )
    log:
        bids(
            root="logs",
            datatype="fit_basic_flatfield",
            subject="{subject}",
            sample="{sample}",
            acq="{acq}",
            stain="{stain}",
            suffix="log.txt",
        ),
    group:
        "preproc"
    container:
        config["containers"]["spimprep"]
    script:
        "../scripts/fit_basic_flatfield_corr_zarr.py"


rule apply_basic_flatfield_corr:
    """ apply BaSiC flatfield correction """
    input:
        zarr=lambda wildcards: bids(
            root=work,
            subject="{subject}",
            datatype="micr",
            sample="{sample}",
            acq="{acq}",
            desc="rawfromgcs" if dataset_is_remote(wildcards) else "raw",
            suffix="SPIM.zarr",
        ).format(**wildcards),
        model_dirs=lambda wildcards: expand(
            rules.fit_basic_flatfield_corr.output.model_dir,
            stain=get_stains(wildcards),
            allow_missing=True,
        ),
    params:
        out_chunks=(1,1,128,256,256) #make this a config option -- setting it here instead of rechunking in zarr2bdv
    output:
        zarr=temp(
            directory(
                bids(
                    root=work,
                    subject="{subject}",
                    datatype="micr",
                    sample="{sample}",
                    acq="{acq}",
                    desc="flatcorr",
                    suffix="SPIM.zarr",
                )
            )
        ),
    benchmark:
        bids(
            root="benchmarks",
            datatype="apply_basic_flatfield",
            subject="{subject}",
            sample="{sample}",
            acq="{acq}",
            suffix="benchmark.tsv",
        )
    log:
        bids(
            root="logs",
            datatype="apply_basic_flatfield",
            subject="{subject}",
            sample="{sample}",
            acq="{acq}",
            suffix="log.txt",
        ),
    resources:
        runtime=60,
        mem_mb=32000,
    threads: config["cores_per_rule"]
    group:
        "preproc"
    script:
        "../scripts/apply_basic_flatfield_corr_zarr.py"
