import json


rule raw_dataset_desc:
    params:
        dd=config["bids"]["raw"],
    output:
        json=Path(root) / "dataset_description.json",
    log:
        "logs/dd_raw.log",
    localrule: True
    run:
        with open(output.json, "w") as outfile:
            json.dump(params.dd, outfile, indent=4)


rule resampled_dataset_desc:
    params:
        dd=config["bids"]["resampled"],
    output:
        json=Path(resampled) / "dataset_description.json",
    log:
        "logs/dd_raw.log",
    localrule: True
    run:
        with open(output.json, "w") as outfile:
            json.dump(params.dd, outfile, indent=4)


rule bids_readme:
    input:
        config["bids"]["readme_md"],
    output:
        Path(root) / "README.md",
    log:
        "logs/bids_readme.log",
    localrule: True
    shell:
        "cp {input} {output}"


rule bids_samples_json:
    input:
        config["bids"]["samples_json"],
    output:
        Path(root) / "samples.json",
    log:
        "logs/bids_samples_json.log",
    localrule: True
    shell:
        "cp {input} {output}"


rule create_samples_tsv:
    input:
        tsv=config["datasets"],
    output:
        tsv=Path(root) / "samples.tsv",
    log:
        "logs/bids_samples_tsv.log",
    localrule: True
    script:
        "../scripts/create_samples_tsv.py"
