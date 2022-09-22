import pandas as pd


boolean_df = pd.read_csv(snakemake.input[0], sep='\t', dtype='object')
scoring_df = pd.read_csv(snakemake.input[1], sep='\t', dtype='object')

not_depleted_col = list(filter(lambda col: col.startswith('Not depleted'), boolean_df.columns))[0]

base_df = boolean_df[[
    'event_id',
    not_depleted_col
]]


features_for_glm = base_df.merge(
    boolean_df[[
        'event_id',
        'algorithms > 1',
        'Type length(mature)',
        'both_donor_acceptor_are_at_annotated_boundary',
        'donor_acceptor_sites_at_the_same_transcript_isoform',
        'both_donor_acceptor_have_AS_events',
    ]],
    on='event_id',
    how='left'
).merge(
    scoring_df[[
        'event_id',
        '#BSJ_reads(totalRNA)',
        'MCS-detail(sample)',
        'evidence_num_plus',
    ]],
    on='event_id',
    how='left'
)[[
    'event_id',
    not_depleted_col,
    '#BSJ_reads(totalRNA)',
    'algorithms > 1',
    'Type length(mature)',
    'MCS-detail(sample)',
    'both_donor_acceptor_are_at_annotated_boundary',
    'donor_acceptor_sites_at_the_same_transcript_isoform',
    'both_donor_acceptor_have_AS_events',
    'evidence_num_plus',
]]


features_for_glm.to_csv(snakemake.output[0], sep='\t', index=False)
