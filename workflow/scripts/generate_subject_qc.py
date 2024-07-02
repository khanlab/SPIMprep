from jinja2 import Environment, FileSystemLoader
import os.path as path

# load jinja template
file_loader = FileSystemLoader(".")
env = Environment(loader=file_loader)
template = env.get_template("qc/resources/subject_html_temp.html")

# input html files
ws_html = snakemake.input.ws_html
ff_html = snakemake.input.ff_html
vol_html = snakemake.input.vol_html

# output html files
sub_html = snakemake.output.sub_html
total_html = "qc/qc_report.html"

# Wildcards
subject = snakemake.wildcards.subject
sample = snakemake.wildcards.sample
acq = snakemake.wildcards.acq

# Get relative path to the subjects QC htmls
ws_rel_path = path.relpath(path.dirname(sub_html), path.dirname(ws_html))+"/"+path.basename(ws_html)
ff_rel_path = path.relpath(path.dirname(sub_html), path.dirname(ff_html))+"/"+path.basename(ff_html)
vol_rel_path = path.relpath(path.dirname(sub_html), path.dirname(vol_html))+ "/" +path.basename(vol_html)

# Fill in jinja template for subject html and write it out
output = template.render(back_link="../qc_report.html",subject=subject,sample=sample,acq=acq,
                         ffhtml=ff_rel_path,wshtml=ws_rel_path, volhtml=vol_rel_path)
with open(sub_html, 'w') as f:
    f.write(output)

# Create line to add link to subject into final qc report combining all subjects
sub_link = f'\n\t\t<a href="./{sub_html.split("/")[1]}/subject.html">{subject}-{sample}-{acq}</a><br>'

# if not first sample just add the one link
if(path.exists(total_html)):
    with open(total_html,'a') as f:
        f.write(sub_link)

# if it is the first sample write out the template
else:
    template = env.get_template("qc/resources/qc_report_temp.html")
    output = template.render()
    output+=sub_link
    with open(total_html, 'w') as f:
        f.write(output)
