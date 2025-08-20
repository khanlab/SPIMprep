#!/usr/bin/env python3
from pathlib import Path

from snakebids import bidsapp, plugins
from plugins import app as app_plugin
from plugins import cli as cli_plugin

app = bidsapp.app(
    [
        app_plugin.SnakemakeBidsApp(Path(__file__).resolve().parent),
        cli_plugin.SpimprepCLIConfig(),
    ]
)


def get_parser():
    """Exposes parser for sphinx doc generation, cwd is the docs dir."""
    return app.build_parser().parser


if __name__ == "__main__":
    app.run()
