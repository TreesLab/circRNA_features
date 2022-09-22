import pandas as pd


circ_df = pd.read_csv(
    snakemake.input[0],
    sep='\t',
    dtype='object',
)

AS_db = pd.read_csv(
    snakemake.input[1],
    sep='\t',
    dtype='object'
)


circ_df_with_AS_check = circ_df.merge(
    AS_db.rename({'pos': 'donor'}, axis=1),
    on=['chr', 'donor', 'strand'],
    how='left'
).merge(
    AS_db.rename({'pos': 'acceptor'}, axis=1),
    on=['chr', 'acceptor', 'strand'],
    how='left',
    suffixes=['(donor)', '(acceptor)']
)

circ_df_with_AS_check.to_csv(snakemake.output[0], sep='\t', index=False)
