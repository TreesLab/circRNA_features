import pandas as pd


circ_df = pd.read_csv(snakemake.input[0], sep='\t', dtype='object')


def append_all_features(circ_df, *feature_dfs):
    for feature_df in feature_dfs:
        circ_df = circ_df.merge(feature_df, on='event_id', how='left')

    return circ_df


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

# donor & acceptor: AS events check
check_AS_df = pd.read_csv(
    snakemake.input.check_AS_events,
    sep='\t',
    dtype='object'
).loc[:, 'event_id':].assign(
    both=lambda df: ((df.iloc[:, 2] == '1') & (df.iloc[:, 4] == '1')).apply(int).apply(str)
).rename(
    {
        'both': 'both_donor_acceptor_have_AS_events'
    },
    axis=1
)

# donor & acceptor: splicing scores
splicing_scores_df = pd.read_csv(
    snakemake.input.splicing_scores,
    sep='\t',
    dtype='object'
).loc[:, 'event_id':]

# donor & acceptor: conservation scores
conservation_scores_df = pd.read_csv(
    snakemake.input.conservation_scores,
    sep='\t',
    dtype='object'
)


# check_ambiguous
check_ambiguous_df = pd.read_csv(
    snakemake.input.check_ambiguous,
    sep='\t',
    dtype='object'
).rename({'circRNA': 'event_id'}, axis=1)


# circAtlas table
circAtlas_df = pd.read_csv(
    snakemake.input.circAtlas,
    sep='\t',
    dtype='object'
)

# transCirc table
transCirc_df = pd.read_csv(
    snakemake.input.transCirc,
    sep='\t',
    dtype='object'
)


# flanking regions: RCS pairs
RCS_df = pd.read_csv(
    snakemake.input.RCS_pairs,
    sep='\t',
    dtype='object'
).rename({'circRNA_id': 'event_id'}, axis=1)

# flanking regions: RBP pairs on flanking 1kb
RBP_df = pd.read_csv(
    snakemake.input.RBP_pairs,
    sep='\t',
    dtype='object',
    names=['event_id', '#RBP_pairs_on_flanking_1k'],
    usecols=[0, 1]
).assign(
    has_common_RBPs_on_flanking_1k=lambda df: (df['#RBP_pairs_on_flanking_1k'].apply(int) > 0).apply(int)
).astype('object')


# circRNA junction: miRNA-binding sites across the circRNA junction
cross_junc_miRNAs_df = pd.read_csv(
    snakemake.input.cross_junc_miRNAs,
    sep='\t',
    dtype='object',
    usecols=[0, 1]
).rename(
    {
        'circRNAs': 'event_id',
        '#miRNAs': '#miRNA_with_binding_sites_across_junction'
    },
    axis=1
).assign(
    has_miRNA_binding_sites_across_junction=lambda df: (df['#miRNA_with_binding_sites_across_junction'].apply(int) > 0).apply(int)
).astype('object')

# circRNA junction: RBP-binding sites across the circRNA junction
cross_junc_RBPs_df = pd.read_csv(
    snakemake.input.cross_junc_RBPs,
    sep='\t',
    dtype='object',
    usecols=[0, 1, 2]
).rename(
    {
        'seq_id': 'event_id',
        'count': '#RBP_with_binding_sites_across_junction',
        'min_p_value': 'min_p_value_of_RBP_binding_sites_across_junction'
    },
    axis=1
).assign(
    has_RBP_binding_sites_across_junction=lambda df: (df['#RBP_with_binding_sites_across_junction'].apply(int) > 0).apply(int)
).astype('object')

cross_junc_RBPs_high_cons_df = pd.read_csv(
    snakemake.input.cross_junc_RBPs_high_cons,
    sep='\t',
    dtype='object',
    usecols=[0, 1, 2]
).rename(
    {
        'seq_id': 'event_id',
        'count': '#RBP_with_binding_sites_across_junction_high_cons',
        'min_p_value': 'min_p_value_of_RBP_binding_sites_across_junction_high_cons'
    },
    axis=1
).assign(
    has_RBP_binding_sites_across_junction_high_cons=lambda df: (df['#RBP_with_binding_sites_across_junction_high_cons'].apply(int) > 0).apply(int)
).astype('object')


# circFL-seq table
circFLseq_df = pd.read_csv(
    snakemake.input.circFLseq,
    sep='\t',
    dtype='object',
    names=['event_id', 'circFLseq_ID']
).assign(
    is_in_circFLseq='1'
).drop('circFLseq_ID', axis=1).drop_duplicates()


# merge all features
circ_df = append_all_features(
    circ_df,
    check_annotated_df,
    check_ambiguous_df,
    check_AS_df,
    circAtlas_df,
    transCirc_df,
    RCS_df,
    RBP_df,
    cross_junc_miRNAs_df,
    cross_junc_RBPs_df,
    cross_junc_RBPs_high_cons_df,
    splicing_scores_df,
    conservation_scores_df,
    circFLseq_df,
).fillna('0')

# output
circ_df.to_csv(snakemake.output[0], sep='\t', index=False)
