import pandas as pd

df = pd.read_csv(snakemake.input[0], sep='\t', dtype='object')

df_with_totalRNA_support = df[df['totalRNA'].apply(int) >= snakemake.params.threshold]
df_with_totalRNA_support.to_csv(snakemake.output[0], sep='\t', index=False)
