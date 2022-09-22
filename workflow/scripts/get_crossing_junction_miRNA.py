import pandas as pd


df = pd.read_csv(snakemake.input[0], sep='\t', dtype='object')
df = df.astype({'ref_start': int, 'ref_end': int, 'score': float})

score_threshold = snakemake.params.score_threshold

df = df[df['score'] >= score_threshold]

L = snakemake.params.L
cross_junction = snakemake.params.cross_junction

cross_df = df[(df['ref_start'] <= L - cross_junction + 1) & (df['ref_end'] >= L + cross_junction)]

miRNA_list = cross_df[[
    'reference_id',
    'query_id'
]].drop_duplicates(
).groupby(
    'reference_id'
).agg(
    [
        'count',
        lambda mirna: ','.join(sorted(mirna))
    ]
).reset_index(
)

miRNA_list.columns = ['circRNAs', '#miRNAs', 'miRNAs']

miRNA_list.to_csv(snakemake.output[0], sep='\t', index=False)
