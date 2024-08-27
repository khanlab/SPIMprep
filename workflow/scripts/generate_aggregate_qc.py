from jinja2 import Environment, FileSystemLoader
from pathlib import Path

datasets = snakemake.params.datasets
total_html = snakemake.output.total_html
subj_htmls = snakemake.input.subj_htmls

# load jinja template
file_loader = FileSystemLoader(".")
env = Environment(loader=file_loader)
template = env.get_template(snakemake.input.report_html)

output = template.render()

for i,subj_html in enumerate(subj_htmls):
    subject=datasets.loc[i,'subject']
    sample=datasets.loc[i,'sample']
    acq=datasets.loc[i,'acq']

    relative_path = Path(subj_html).relative_to(Path(total_html).parent)
    # Create line to add link to subject into final qc report combining all subjects
    subj_link = f'\n\t\t<a href="{relative_path}">sub-{subject}_sample-{sample}_acq-{acq}</a><br>\n'

    output+=subj_link

with open(total_html, 'w') as f:
    f.write(output)
