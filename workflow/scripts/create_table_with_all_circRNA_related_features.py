import pandas as pd


circ_df = pd.read_csv(snakemake.input.circRNAs, sep='\t', dtype='object')
event_id_df = circ_df[['event_id']]

def append_all_features(circ_df, *feature_dfs):
    for feature_df in feature_dfs:
        circ_df = circ_df.merge(feature_df, on='event_id', how='left')

    return circ_df



# not depleted ratio
not_depleted_ratio_df = pd.read_csv(
    snakemake.input.not_depleted_ratio,
    sep='\t',
    dtype='object'
).merge(event_id_df, on='event_id', how='right')


# circFL-seq table
circFLseq_df = pd.read_csv(
    snakemake.input.circFLseq,
    sep='\t',
    dtype='object',
    names=['event_id', 'circFLseq_ID']
).assign(
    is_in_circFLseq='1'
).drop(
    'circFLseq_ID',
    axis=1
).drop_duplicates().merge(event_id_df, on='event_id', how='outer').fillna('0')


# circAtlas table: Type Length, Length
circAtlas_df = pd.read_csv(
    snakemake.input.circAtlas,
    sep='\t',
    dtype='object'
)[[
    'event_id',
    'Type length',
    'Length (nt)'
]]


# donor & acceptor: annotation check
check_annotated_df = pd.read_csv(
    snakemake.input.check_annotated,
    sep='\t',
    dtype='object',
    usecols=[4, 5, 6, 7, 8]
).rename(
    {
        'circ_id': 'event_id',
        'host_gene': 'host_gene(circmimi)'
    },
    axis=1
).assign(
    both=lambda df: ((df.iloc[:, 2] == '1') & (df.iloc[:, 3] == '1')).apply(int).apply(str)
).rename(
    {
        'both': 'both_donor_acceptor_are_at_annotated_boundary'
    },
    axis=1
)


# check_ambiguous
check_ambiguous_df = pd.read_csv(
    snakemake.input.check_ambiguous,
    sep='\t',
    dtype='object'
).rename({'circRNA': 'event_id'}, axis=1)


# merge all features
circ_df = append_all_features(
    circ_df,
    circAtlas_df,
    circFLseq_df,
    check_annotated_df,
    check_ambiguous_df,
    not_depleted_ratio_df,
)

# output
circ_df.to_csv(snakemake.output[0], sep='\t', index=False)
