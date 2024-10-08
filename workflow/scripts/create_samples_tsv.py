import pandas as pd

df = snakemake.params.samples_df

df['participant_id'] = 'sub-' + df['subject'] 
df['sample_id'] = 'sample-' + df['sample'] 
df['sample_type'] = 'tissue'

df[['participant_id','sample_id','sample_type']].to_csv(snakemake.output.tsv,sep='\t',index=False)


