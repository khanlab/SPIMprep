import csv
import subprocess
from pathlib import Path
import click
from rich.console import Console
from rich.table import Table
from rich.prompt import Confirm
from pathlib import Path

# Path to this script
SCRIPT_DIR = Path(__file__).resolve().parent

# Path to run.py (assumed to be in the same dir)
RUN_PY_PATH = SCRIPT_DIR / "run.py"


console = Console()


def submit_job(row, sbatch_script, log_dir, slurm_config,
               output_bids_dir, work_dir, conda_prefix, dry_run=False):

    subject = row["subject"]
    acq = row["acq"]
    stain_0 = row['stain_0']
    stain_1 = row['stain_1']
    stain_2 = row['stain_2']
    input_path = row["sample_path"]

    job_name = f"{subject}_{acq}"

    log_out = log_dir / f"{job_name}.out"

    slurm_args = [
        "--job-name", job_name,
        "--time", slurm_config["time"],
        "--mem", slurm_config["mem"],
        "--nodelist", "rri-cbs-h2.schulich.uwo.ca",
        "--cpus-per-task", slurm_config["cpus"],
        "--tmp", slurm_config["tmp"],
        "--output", str(log_out),
        "--export", f"ALL,SUBJECT={subject},ACQ={acq},STAIN_0={stain_0},STAIN_1={stain_1},STAIN_2={stain_2},INPUT_PATH={input_path},OUTPUT_BIDS_DIR={output_bids_dir},WORK_DIR={work_dir},CONDA_PREFIX={conda_prefix},RUN_PY_PATH={RUN_PY_PATH}"
    ]

    cmd = ["sbatch"] + slurm_args + [str(sbatch_script)]

    if dry_run:
        console.rule(f"[bold yellow]Dry run: {job_name}[/bold yellow]")

        console.print("[cyan]sbatch command:[/cyan]")
        console.print(" ".join(cmd))

        console.print("\n[magenta]Exported environment variables:[/magenta]")
        for var in ["SUBJECT", "ACQ", "STAIN_0", "STAIN_1", "STAIN_2", "INPUT_PATH", "OUTPUT_BIDS_DIR", "WORK_DIR", "CONDA_PREFIX"]:
            value = eval(var.lower())  # grab from local variables
            console.print(f"{var} = \"{value}\"")

        console.print("\n[green]Simulated input line in run_sample.sh:[/green]")
        console.print(f'--input-path "{input_path}"')  # confirm proper quoting

        console.rule()
    else:
        console.print(f"[green]Submitting:[/green] [bold]{job_name}[/bold]")
        subprocess.run(cmd)


@click.command()
@click.option('--samples-file', '-s', type=click.Path(exists=True,dir_okay=False, file_okay=True, readable=True), default="config/samples.tsv", 
              help="Path to the TSV file containing sample metadata.")
@click.option('--sbatch-script', '-b', type=click.Path(exists=True,dir_okay=False, file_okay=True, readable=True), default="run_sample.sh",
              help="Path to the SBATCH script template.")
@click.option('--log-dir', '-l', type=click.Path(), default="logs",
              help="Directory where SLURM output logs are stored.")
@click.option('--output-bids-dir', required=True, type=click.Path(), help="Path to BIDS output directory.")
@click.option('--work-dir', required=True, type=click.Path(exists=True,dir_okay=True, file_okay=False, readable=True), help="Path to scratch or work directory.", default="/tmp")
@click.option('--conda-prefix', required=True, type=click.Path(exists=True,dir_okay=True, file_okay=False, readable=True), help="Path to conda prefix for Snakemake.")
@click.option('--time', default="2-00:00:00", help="SLURM time allocation per job.")
@click.option('--mem', default="64G", help="SLURM memory per job.")
@click.option('--cpus', default="16", help="SLURM CPUs per job.")
@click.option('--tmp', default="4T", help="Temporary disk space per job.")
@click.option('--dry-run', is_flag=True, help="Only print the sbatch commands without running.")
@click.option('--yes', '-y', is_flag=True, help="Skip confirmation before submitting jobs.")
def main(samples_file, sbatch_script, log_dir,
         output_bids_dir, work_dir, conda_prefix,
         time, mem, cpus, tmp, dry_run, yes):
    """Submit SPIMprep jobs for each sample in a TSV file via sbatch."""

    log_dir = Path(log_dir)
    log_dir.mkdir(exist_ok=True)

    slurm_config = {"time": time, "mem": mem, "cpus": cpus, "tmp": tmp}

    with open(samples_file, newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        rows = list(reader)

    if not yes:
        table = Table(title="Jobs to Submit", show_lines=True)
        table.add_column("Subject")
        table.add_column("Acq")
        table.add_column("Stains")
        table.add_column("Input Path")

        for row in rows:
            stains = f"{row['stain_0']} {row['stain_1']} {row['stain_2']}"
            table.add_row(row["subject"], row["acq"], stains, row["sample_path"])

        console.print(table)
        if not Confirm.ask("Submit these jobs?"):
            console.print("[red]Cancelled.[/red]")
            return

    for row in rows:
        submit_job(
            row=row,
            sbatch_script=Path(sbatch_script),
            log_dir=log_dir,
            slurm_config=slurm_config,
            dry_run=dry_run,
            output_bids_dir=output_bids_dir,
            work_dir=work_dir,
            conda_prefix=conda_prefix
        )


if __name__ == "__main__":
    main()

