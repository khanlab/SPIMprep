
rule setup_qc_dir:
    "Copies QC resources to the output bids folder"
    input:
        readme_md=config["report"]["resources"]["readme_md"],
    output:
        readme_md=remote_file(Path(root) / "qc" / "README.md"),
    threads: 1
    resources:
        mem_mb=1000,
        runtime=10,
    log:
        "logs/setup_qc_dir_log.txt",
    shell:
        "cp {input.readme_md} {output.readme_md}"


rule generate_flatfield_qc:
    "Generates an html file for comparing before and after flatfield correction"
    input:
        uncorr=lambda wildcards: bids(
            root=work,
            subject="{subject}",
            datatype="micr",
            sample="{sample}",
            acq="{acq}",
            desc="rawfromgcs" if sample_is_remote(wildcards) else "raw",
            suffix="SPIM.zarr",
        ).format(**wildcards),
        corr=bids(
            root=work,
            subject="{subject}",
            datatype="micr",
            sample="{sample}",
            acq="{acq}",
            desc="flatcorr",
            suffix="SPIM.zarr",
        ),
        ff_html=config["report"]["resources"]["ff_html"],
    params:
        ff_s_start=config["report"]["flatfield_corrected"]["slice_start"],
        ff_s_step=config["report"]["flatfield_corrected"]["slice_step"],
        ff_cmap=config["report"]["flatfield_corrected"]["colour_map"],
    output:
        html=remote_file(
            Path(root)
            / "qc"
            / "sub-{subject}_sample-{sample}_acq-{acq}/flatfieldqc.html"
        ),
        corr_images_dir=remote_directory(
            Path(root)
            / "qc"
            / "sub-{subject}_sample-{sample}_acq-{acq}"
            / "images"
            / "corr"
        ),
        uncorr_images_dir=remote_directory(
            Path(root)
            / "qc"
            / "sub-{subject}_sample-{sample}_acq-{acq}"
            / "images"
            / "uncorr"
        ),
    threads: 1
    resources:
        mem_mb=8000,
        runtime=60,
    log:
        bids(
            root="logs",
            datatype="generate_flatfield_qc",
            subject="{subject}",
            sample="{sample}",
            acq="{acq}",
            suffix="log.txt",
        ),
    script:
        "../scripts/generate_flatfield_qc.py"


rule generate_whole_slice_qc:
    "Generates an html file to view whole slices from preprocessed images"
    input:
        **get_storage_creds(),
        zarr=bids(
            root=root,
            subject="{subject}",
            datatype="micr",
            sample="{sample}",
            acq="{acq}",
            suffix="SPIM.{ext}".format(ext=get_extension_ome_zarr()),
        ),
        ws_html=config["report"]["resources"]["ws_html"],
    params:
        ws_s_start=config["report"]["whole_slice_viewer"]["slice_start"],
        ws_s_step=config["report"]["whole_slice_viewer"]["slice_step"],
        ws_cmap=config["report"]["whole_slice_viewer"]["colour_map"],
        uri=get_output_ome_zarr_uri(),
        storage_provider_settings=workflow.storage_provider_settings,
    output:
        html=remote_file(
            Path(root)
            / "qc"
            / "sub-{subject}_sample-{sample}_acq-{acq}"
            / "whole_slice_qc.html"
        ),
        images_dir=remote_directory(
            Path(root)
            / "qc"
            / "sub-{subject}_sample-{sample}_acq-{acq}"
            / "images"
            / "whole"
        ),
    threads: 1
    resources:
        mem_mb=8000,
        runtime=60,
    log:
        bids(
            root="logs",
            datatype="generate_whole_slice_qc",
            subject="{subject}",
            sample="{sample}",
            acq="{acq}",
            suffix="log.txt",
        ),
    script:
        "../scripts/generate_whole_slice_qc.py"


rule generate_volume_qc:
    "Generates an html file to view the volume rendered image"
    input:
        **get_storage_creds(),
        zarr=bids(
            root=root,
            subject="{subject}",
            datatype="micr",
            sample="{sample}",
            acq="{acq}",
            suffix="SPIM.{ext}".format(ext=get_extension_ome_zarr()),
        ),
        vol_viewer_dir=config["report"]["resources"]["vol_viewer_dir"],
    params:
        uri=get_output_ome_zarr_uri(),
        storage_provider_settings=workflow.storage_provider_settings,
    output:
        resources=remote_directory(
            Path(root)
            / "qc"
            / "sub-{subject}_sample-{sample}_acq-{acq}"
            / "volume_resources"
        ),
        html=remote_file(
            Path(root)
            / "qc"
            / "sub-{subject}_sample-{sample}_acq-{acq}"
            / "volume_qc.html"
        ),
    threads: 1
    resources:
        mem_mb=8000,
        runtime=60,
    log:
        bids(
            root="logs",
            datatype="generate_volume_qc",
            subject="{subject}",
            sample="{sample}",
            acq="{acq}",
            suffix="log.txt",
        ),
    script:
        "../scripts/generate_volume_qc.py"


rule generate_subject_qc:
    "Generates html files to access all the subjects qc reports in one place"
    input:
        subject_html=config["report"]["resources"]["subject_html"],
        ws_html=rules.generate_whole_slice_qc.output.html,
        ff_html=rules.generate_flatfield_qc.output.html,
        vol_html=rules.generate_volume_qc.output.html,
    output:
        sub_html=remote_file(
            Path(root)
            / "qc"
            / "sub-{subject}_sample-{sample}_acq-{acq}"
            / "subject.html"
        ),
    threads: 1
    resources:
        mem_mb=2000,
        runtime=10,
    log:
        bids(
            root="logs",
            datatype="generate_subject_qc",
            subject="{subject}",
            sample="{sample}",
            acq="{acq}",
            suffix="log.txt",
        ),
    script:
        "../scripts/generate_subject_qc.py"


rule generate_aggregate_qc:
    input:
        report_html=config["report"]["resources"]["report_html"],
        subj_htmls=get_all_subj_html,
    params:
        samples=samples,
    output:
        total_html=remote_file(Path(root) / "qc" / "qc_report.html"),
    threads: 1
    resources:
        mem_mb=2000,
        runtime=10,
    log:
        bids(
            root="logs",
            datatype="generate_aggregate_qc",
            suffix="log.txt",
        ),
    script:
        "../scripts/generate_aggregate_qc.py"
