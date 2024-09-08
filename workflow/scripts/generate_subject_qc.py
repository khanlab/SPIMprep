from jinja2 import Environment, FileSystemLoader
import os.path as path
from pathlib import Path

# load jinja template
file_loader = FileSystemLoader(".")
env = Environment(loader=file_loader)
template = env.get_template(snakemake.input.subject_html)

# input html files
in_htmls = snakemake.input.in_htmls
#ws_html = snakemake.input.ws_html
#ff_html = snakemake.input.ff_html
#vol_html = snakemake.input.vol_html

# output html files
sub_html = snakemake.output.sub_html

# Wildcards
subject = snakemake.wildcards.subject
sample = snakemake.wildcards.sample
acq = snakemake.wildcards.acq

# Get relative path to the subjects QC htmls

link_rel_paths = [ Path(link_path).relative_to(Path(sub_html).parent) for link_path in snakemake.input.link_paths]
link_descriptions = snakemake.params.link_descriptions

# Fill in jinja template for subject html and write it out
output = template.render(back_link="../qc_report.html",subject=subject,sample=sample,acq=acq,
                         link_rel_paths, link_descriptions)
with open(sub_html, 'w') as f:
    f.write(output)

