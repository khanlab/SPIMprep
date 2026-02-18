from __future__ import annotations

import argparse
from pathlib import Path
from tempfile import gettempdir
from typing import Any

import attrs
from snakebids import bidsapp
from snakebids.bidsapp.args import ArgumentGroups
from snakebids.plugins.base import PluginBase


@attrs.define
class SpimprepCLIConfig(PluginBase):
    """Dynamically add CLI parameters for SPIMprep.

    Parameters
    ----------
    argument_group
        Specify title of the group to which arguments should be added


    CLI Arguments
    ~~~~~~~~~~~~~
    Adds a atlas argument to the CLI, along with the related configuration
    for the atlas into the config.
    """

    argument_group: str | None = None

    @bidsapp.hookimpl
    def add_cli_arguments(
        self, parser: argparse.ArgumentParser, argument_groups: ArgumentGroups
    ):
        """Add cli parameters."""
        group = (
            argument_groups[self.argument_group]
            if self.argument_group is not None
            else parser
        )
        self.try_add_argument(
            group,
            "--work-dir",
            action="store",
            type=str,
            dest="work_dir",
            default=None,
            help=(
                "Set the work directory (effectively the snakebids output workflow directory)"
            ),
        )
        self.try_add_argument(
            group,
            "--output-bids-dir",
            action="store",
            type=str,
            dest="output_bids_dir",
            required=True,
            help=("Set the output bids directory"),
        )

        self.try_add_argument(
            group,
            "--stains",
            action="store",
            type=str,
            nargs="+",
            dest="stains",
            required=True,
            help=("Set the stains for each channels"),
        )
        self.try_add_argument(
            group,
            "--subject",
            action="store",
            type=str,
            dest="subject",
            required=True,
            help=("Set the subject identifier (participant-label)"),
        )
        self.try_add_argument(
            group,
            "--acq",
            action="store",
            type=str,
            default="blaze",
            dest="acq",
            help=("Set the acquisition entity"),
        )
        self.try_add_argument(
            group,
            "--sample",
            action="store",
            type=str,
            default="brain",
            dest="sample",
            help=("Set the sample entity"),
        )
        self.try_add_argument(
            group,
            "--input-path",
            action="store",
            type=str,
            dest="input_path",
            required=True,
            help=("Set the input path"),
        )
        self.try_add_argument(
            group,
            "--resolution-preset",
            action="store",
            type=str,
            dest="resolution_preset",
            default=None,
            help=(
                "Use a named resolution preset from config to override "
                "physical sizes and tif glob for prestitched imports. "
                "Available presets are defined in config.yml under "
                "import_prestitched.resolution_presets"
            ),
        )

    @bidsapp.hookimpl
    def update_cli_namespace(self, namespace: dict[str, Any], config: dict[str, Any]):
        """Add vars to config."""

        work_dir = self.pop(namespace, "work_dir")

        if work_dir == None:
            work_dir = gettempdir()

        config["work_dir"] = Path(work_dir)

        output_bids_dir = self.pop(namespace, "output_bids_dir")

        config["output_dir"] = Path(output_bids_dir)

        # Apply resolution preset if specified
        resolution_preset = self.pop(namespace, "resolution_preset")
        if resolution_preset is not None:
            presets = config.get("import_prestitched", {}).get(
                "resolution_presets", {}
            )
            if resolution_preset not in presets:
                available = ", ".join(presets.keys()) if presets else "(none defined)"
                raise ValueError(
                    f"Unknown resolution preset '{resolution_preset}'. "
                    f"Available presets: {available}"
                )
            preset = presets[resolution_preset]
            for key, value in preset.items():
                config["import_prestitched"][key] = value

    @bidsapp.hookimpl
    def finalize_config(self, config: dict[str, Any]) -> None:
        """create the samples.tsv required for the workflow"""

        # Prepare samples.tsv
        sample_info = {
            "subject": config["subject"],
            "sample": config["sample"],
            "acq": config["acq"],
            "stain_0": config["stains"][0],
            "sample_path": config["input_path"],
        }
        for i, stain_var in enumerate(config["stains"][1:], start=1):
            sample_info[f"stain_{i}"] = stain_var

        samples_tsv_path = config["work_dir"] / "samples_from_cli.tsv"
        with open(samples_tsv_path, "w") as f:
            headers = "\t".join(sample_info.keys())
            f.write(headers + "\n")
            values = "\t".join(sample_info.values())
            f.write(values + "\n")

        config["samples"] = "samples_from_cli.tsv"