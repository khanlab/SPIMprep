rule raw_dataset_desc:
    params:
        dd=config["bids"]["raw"],
    output:
        json=Path(root) / "dataset_description.json",
    log:
        "logs/dd_raw.log",
    run:
        import json

        with open(output.json, "w") as outfile:
            json.dump(params.dd, outfile, indent=4)


rule resampled_dataset_desc:
    params:
        dd=config["bids"]["resampled"],
    output:
        json=Path(resampled) / "dataset_description.json",
    log:
        "logs/dd_raw.log",
    run:
        import json

        with open(output.json, "w") as outfile:
            json.dump(params.dd, outfile, indent=4)
