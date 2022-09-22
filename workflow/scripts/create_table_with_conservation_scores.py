import pandas as pd


df = pd.read_csv(
    snakemake.input[0],
    sep='\t',
    dtype='object',
    names=['chr', 'pos1', 'pos2', 'event_id', 'type', 'phyloP', 'phastCons']
).loc[:, 'event_id':]

ev_id_df = df[['event_id']].drop_duplicates().reset_index(drop=True)

unstack_df = df.set_index(['event_id', 'type']).unstack(level=1).reset_index()
unstack_df.columns = ['event_id'] + ["{}({})".format(*col) for col in unstack_df.columns[1:]]

unstack_df.to_csv(snakemake.output[0],  sep='\t', index=False)
