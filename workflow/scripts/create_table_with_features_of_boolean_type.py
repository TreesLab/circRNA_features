import pandas as pd
import logging


logging.basicConfig(filename=snakemake.log[0], encoding='utf-8', level=logging.DEBUG)


logging.info('loading raw features table')

df = pd.read_csv(
    snakemake.input[0],
    sep='\t',
    dtype='object'
).loc[:, 'event_id':].set_index('event_id')


# not detected & not depleted
logging.info('generating binary feature of "not detected & not depleted" columns')

RNaseR_polyAminus_col = list(filter(lambda col: col in ['RNaseR', 'polyA_minus'], snakemake.params.treatments))[0]
logging.debug(RNaseR_polyAminus_col)
not_detected_col = f'Not detected({RNaseR_polyAminus_col})'
not_depleted_col = f'Not depleted({RNaseR_polyAminus_col})'
logging.debug([not_detected_col, not_depleted_col])
logging.debug(df.columns)

base_df = df[[not_detected_col, not_depleted_col]]


# donor & acceptor: annotation check
logging.info('generating binary feature of "donor & acceptor: annotation check"')

annotated_boundary_df = df[[
    'donor_site_at_the_annotated_boundary',
    'acceptor_site_at_the_annotated_boundary',
    'both_donor_acceptor_are_at_annotated_boundary',
    'donor_acceptor_sites_at_the_same_transcript_isoform',
]]

# donor & acceptor: AS events check
logging.info('generating binary feature of "donor & acceptor: AS events check"')

AS_check_df = df[[
    'has_AS_event(donor)',
    'has_AS_event(acceptor)',
    'both_donor_acceptor_have_AS_events',
]]


# RCS
logging.info('generating binary feature of "RCS-related features"')

RCS_df = df[[
    '#RCS_across_flanking_introns',
    '#RCS_within_donor\'s_intron',
    '#RCS_within_acceptor\'s_intron'
]].fillna('0').astype(int).assign(
    **{
        '#RCS_across > 0': lambda df2: (df2['#RCS_across_flanking_introns'] > 0).apply(int),
        '#RCS_across - #RCS_within > 0': lambda df2: ((df2['#RCS_across_flanking_introns'] - df2['#RCS_within_donor\'s_intron'] - df2['#RCS_within_acceptor\'s_intron']) > 0).apply(int)
    }
).drop(
    [
        '#RCS_across_flanking_introns',
        '#RCS_within_donor\'s_intron',
        '#RCS_within_acceptor\'s_intron'
    ],
    axis=1
)

# circAtlas table: Type length(mature)
logging.info('generating binary feature of "circAtlas table: Type length"')

type_length_df = df[[
    'Type length(mature)',
    'Type length(Estimated)',
    'Type length(Unknown)'
]]


# circAtlas table: MCS
logging.info('generating binary feature of "MCS"')

MCS_df = df[[
    'MCS-detail(species)',
    'MCS-detail(tissue)',
    'MCS-detail(sample)'
]].astype(int).assign(
    **{
        'MCS-detail(species) > 1': lambda df2: (df2['MCS-detail(species)'] > 1).apply(int),
        'MCS-detail(species) > 2': lambda df2: (df2['MCS-detail(species)'] > 2).apply(int),
        'MCS-detail(species) > 3': lambda df2: (df2['MCS-detail(species)'] > 3).apply(int),
        'MCS-detail(tissue) > 1': lambda df2: (df2['MCS-detail(tissue)'] > 1).apply(int),
        'MCS-detail(tissue) > 2': lambda df2: (df2['MCS-detail(tissue)'] > 2).apply(int),
        'MCS-detail(tissue) > 3': lambda df2: (df2['MCS-detail(tissue)'] > 3).apply(int),
        'MCS-detail(sample) > 1': lambda df2: (df2['MCS-detail(sample)'] > 1).apply(int),
        'MCS-detail(sample) > 2': lambda df2: (df2['MCS-detail(sample)'] > 2).apply(int),
        'MCS-detail(sample) > 3': lambda df2: (df2['MCS-detail(sample)'] > 3).apply(int)
    }
).drop(
    [
        'MCS-detail(species)',
        'MCS-detail(tissue)',
        'MCS-detail(sample)'
    ],
    axis=1
)

# circAtlas table: algorithms
logging.info('generating binary feature of "algorithms"')

algorithms_df = df[[
    'Algorithm(CIRI2)',
    'Algorithm(DCC)',
    'Algorithm(find_circ)',
    'Algorithm(CIRCexplorer2)'
]].astype(int).assign(
    algorithms=lambda df2: df2.apply(sum, axis=1)
).assign(
    **{
        'algorithms > 1': lambda df2: (df2['algorithms'] > 1).apply(int),
        'algorithms > 2': lambda df2: (df2['algorithms'] > 2).apply(int),
        'algorithms > 3': lambda df2: (df2['algorithms'] > 3).apply(int),
    }
).drop(
    [
        'Algorithm(CIRI2)',
        'Algorithm(DCC)',
        'Algorithm(find_circ)',
        'Algorithm(CIRCexplorer2)',
        'algorithms'
    ],
    axis=1
)


# transCirc table: evidence_num > N
logging.info('generating binary feature of "transCirc table: evidence_num > N"')

evidence_num_df = df[[
    'evidences_num'
]].fillna('0').astype(int).assign(
    **{
        'evidences_num > 1': lambda df2: (df2['evidences_num'] > 1).apply(int),
        'evidences_num > 2': lambda df2: (df2['evidences_num'] > 2).apply(int),
        'evidences_num > 3': lambda df2: (df2['evidences_num'] > 3).apply(int),
    }
).drop('evidences_num', axis=1)


# others
logging.info('generating binary feature of the others')

the_others_df = df[[
    'has_common_RBPs_on_flanking_1k',
    'has_miRNA_binding_sites_across_junction',
    'has_RBP_binding_sites_across_junction',
    'has_RBP_binding_sites_across_junction_high_cons',
    'is_in_circFLseq'
]]


# G quadruplex structure across circRNA junction
has_G_quadruplex_10_df = df[[
    'G score(d=10)'
]].astype(int).assign(
    has_G_quadruplex_across_junction_10=lambda sdf: (sdf['G score(d=10)'] > 0).apply(int)
).drop(
    'G score(d=10)',
    axis=1
)

has_G_quadruplex_5_df = df[[
    'G score(d=5)'
]].astype(int).assign(
    has_G_quadruplex_across_junction_5=lambda sdf: (sdf['G score(d=5)'] > 0).apply(int)
).drop(
    'G score(d=5)',
    axis=1
)


logging.info('merging all feature tables')

bool_table_df = pd.concat(
    [
        base_df,
        annotated_boundary_df,
        AS_check_df,
        RCS_df,
        type_length_df,
        MCS_df,
        algorithms_df,
        evidence_num_df,
        the_others_df,
        has_G_quadruplex_10_df,
        has_G_quadruplex_5_df
    ],
    axis=1
)

logging.info('saving results...')

bool_table_df.reset_index().to_csv(snakemake.output[0], sep='\t', index=False)

logging.info('done!')
