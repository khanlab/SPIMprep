rule generate_flatfield_qc:
    "Generates an html file for comparing before and after flatfield correction"
    input:
        uncorr=bids(
            root=work,
            subject="{subject}",
            datatype="micr",
            sample="{sample}",
            acq="{acq}",
            desc="raw",
            suffix="SPIM.zarr",
        ),
        corr=bids(
            root=work,
            subject="{subject}",
            datatype="micr",
            sample="{sample}",
            acq="{acq}",
            desc="flatcorr",
            suffix="SPIM.zarr",
        ),
    params:
        ff_s_start=config["report"]["flatfield_corrected"]["slice_start"],
        ff_s_step=config["report"]["flatfield_corrected"]["slice_step"],
        ff_cmap=config["report"]["flatfield_corrected"]["colour_map"],
    output:
        html="qc/sub-{subject}_sample-{sample}_acq-{acq}/flatfieldqc.html",
        corr_images_dir=directory(
            "qc/sub-{subject}_sample-{sample}_acq-{acq}/images/corr"
        ),
        uncorr_images_dir=directory(
            "qc/sub-{subject}_sample-{sample}_acq-{acq}/images/uncorr"
        ),
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
        ome=get_input_ome_zarr_to_nii(),
    params:
        ws_s_start=config["report"]["whole_slice_viewer"]["slice_start"],
        ws_s_step=config["report"]["whole_slice_viewer"]["slice_step"],
        ws_cmap=config["report"]["whole_slice_viewer"]["colour_map"],
    output:
        html="qc/sub-{subject}_sample-{sample}_acq-{acq}/whole_slice_qc.html",
        images_dir=directory("qc/sub-{subject}_sample-{sample}_acq-{acq}/images/whole"),
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
        ome=get_input_ome_zarr_to_nii(),
    output:
        resources=directory(
            "qc/sub-{subject}_sample-{sample}_acq-{acq}/volume_resources"
        ),
        html="qc/sub-{subject}_sample-{sample}_acq-{acq}/volume_qc.html",
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
        ws_html="qc/sub-{subject}_sample-{sample}_acq-{acq}/whole_slice_qc.html",
        ff_html="qc/sub-{subject}_sample-{sample}_acq-{acq}/flatfieldqc.html",
        vol_html="qc/sub-{subject}_sample-{sample}_acq-{acq}/volume_qc.html",
    output:
        sub_html="qc/sub-{subject}_sample-{sample}_acq-{acq}/subject.html",
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
