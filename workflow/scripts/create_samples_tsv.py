import pandas as pd

# read datasets tsv
df = pd.read_csv(
    snakemake.input.tsv,
    sep="\t",
    dtype={
        "subject": str,
        "sample": str,
        "acq": str,
        "stain_0": str,
        "stain_1": str,
        "num_tiles": int,
    },
)

df['participant_id'] = 'sub-' + df['subject'] 
df['sample_id'] = 'sample-' + df['sample'] 
df['sample_type'] = 'tissue'

df[['participant_id','sample_id','sample_type']].to_csv(snakemake.output.tsv,sep='\t',index=False)


