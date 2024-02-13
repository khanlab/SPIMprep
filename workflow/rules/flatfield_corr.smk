
rule fit_basic_flatfield_corr:
    """ BaSiC flatfield correction"""
    input:
        zarr=bids(
            root=work,
            subject="{subject}",
            datatype="micr",
            sample="{sample}",
            acq="{acq}",
            desc="raw",
            suffix="spim.zarr",
        ),
    params:
        channel=lambda wildcards: get_stains(wildcards).index(wildcards.stain),
        max_n_images=config["basic_flatfield_corr"]["max_n_images"],  #sets maximum number of images to use for fitting (selected randomly)
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
        runtime=90,  #this should be proportional to the number of images and image size
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
    group:
        "preproc"
    container:
        config["containers"]["spimprep"]
    script:
        "../scripts/fit_basic_flatfield_corr_zarr.py"


rule apply_basic_flatfield_corr:
    """ apply BaSiC flatfield correction """
    input:
        zarr=bids(
            root=work,
            subject="{subject}",
            datatype="micr",
            sample="{sample}",
            acq="{acq}",
            desc="raw",
            suffix="spim.zarr",
        ),
        model_dirs=lambda wildcards: expand(
            rules.fit_basic_flatfield_corr.output.model_dir,
            stain=get_stains(wildcards),
            allow_missing=True,
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
                    desc="flatcorr",
                    suffix="spim.zarr",
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
    resources:
        runtime=60,  #this should be proportional to the number of images and image size
        mem_mb=32000,
    threads: 32
    group:
        "preproc"
    script:
        "../scripts/apply_basic_flatfield_corr_zarr.py"
