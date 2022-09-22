import pandas as pd

df = pd.read_csv(snakemake.input[0], sep='\t', dtype='object')

no_AA_df = df[
    (df['totalRNA'].astype(int) >= int(snakemake.wildcards.totalRNA_threshold)) & 
    (df['colinear'] == '0') & 
    (df['multiple_hits'] == '0')
]
no_AA_df.to_csv(snakemake.output[0], sep='\t', index=False)
