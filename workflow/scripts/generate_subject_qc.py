import os.path as path
from pathlib import Path

from jinja2 import Environment, FileSystemLoader

# load jinja template
file_loader = FileSystemLoader(".")
env = Environment(loader=file_loader)
template = env.get_template(snakemake.input.subject_html)

# input html files
ws_html = snakemake.input.ws_html
ff_html = snakemake.input.ff_html
vol_html = snakemake.input.vol_html

# output html files
sub_html = snakemake.output.sub_html

# Wildcards
subject = snakemake.wildcards.subject
sample = snakemake.wildcards.sample
acq = snakemake.wildcards.acq

# Get relative path to the subjects QC htmls

ws_rel_path = Path(ws_html).relative_to(Path(sub_html).parent)
ff_rel_path = Path(ff_html).relative_to(Path(sub_html).parent)
vol_rel_path = Path(vol_html).relative_to(Path(sub_html).parent)

# Fill in jinja template for subject html and write it out
output = template.render(
    back_link="../qc_report.html",
    subject=subject,
    sample=sample,
    acq=acq,
    ffhtml=ff_rel_path,
    wshtml=ws_rel_path,
    volhtml=vol_rel_path,
)
with open(sub_html, "w") as f:
    f.write(output)
