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
not_detected_col = f'Not detected({RNaseR_polyAminus_col})'
not_depleted_col = f'Not depleted({RNaseR_polyAminus_col})'

base_df = df[[not_detected_col, not_depleted_col]]


# numReads(totalRNA)
num_reads_df = df[[
    'totalRNA'
]].rename(
    {
        'totalRNA': '#BSJ_reads(totalRNA)'
    },
    axis=1
)

# RPM
RPM_totalRNA_col = list(filter(lambda col: col.startswith('RPM_totalRNA'), df.columns))[0]
RPM_df = df[[RPM_totalRNA_col]].rename({RPM_totalRNA_col: 'RPM_totalRNA'}, axis=1)


# circAtlas table: MCS
logging.info('generating scoring feature of "circAtlas table: MCS"')

MCS_df = df[[
    'Multiple Conservation Score(MCS)',
    'MCS-detail(species)',
    'MCS-detail(tissue)',
    'MCS-detail(sample)'
]]

# circAtlas table: Tissue Specifity Index
logging.info('generating scoring feature of "circAtlas table: Tissue Specifity Index"')

TSI_df = df[[
    'Tissue Specifity Index'
]]


# transCirc table: evidence
logging.info('generating scoring feature of "transCirc table: evidence"')

evidence_df = df[[
    'evidences_num',
    'evidences_score'
]].fillna('0').astype(
    {
        'evidences_num': int,
        'evidences_score': float
    }
)

# evidence_num + "miRNA/RBP cross junction" 
evidence_num_plus_df = df[[
    'evidences_num',
    'has_miRNA_binding_sites_across_junction',
    'has_RBP_binding_sites_across_junction_high_cons',
]].fillna('0').astype(int).assign(
    evidence_num_plus=lambda sdf: sdf.sum(axis=1)
)[[
    'evidence_num_plus'
]]


# flanking_regions: 

flanking_regions_df = df[[
    '#RCS_across_flanking_introns',
    'min_dist_of_RCS(across, donor)',
    'min_dist_of_RCS(across, acceptor)',
    'sum_of_the_min_dist',
    '#RCS_within_donor\'s_intron',
    '#RCS_within_acceptor\'s_intron',
    '#RBP_pairs_on_flanking_1k',
    '#miRNA_with_binding_sites_across_junction',
    '#RBP_with_binding_sites_across_junction',
    '#RBP_with_binding_sites_across_junction_high_cons'
]].fillna('0').astype(int).assign(
    **{
        '#RCS_across - #RCS_within': lambda df2: (df2['#RCS_across_flanking_introns'] - df2['#RCS_within_donor\'s_intron'] - df2['#RCS_within_acceptor\'s_intron'])
    }
)

flanking_regions_df_2 = df[[
    'min_p_value_of_RBP_binding_sites_across_junction',
    'min_p_value_of_RBP_binding_sites_across_junction_high_cons'
]].fillna('0').astype(float)


# splicing scores

splicing_scores_df = df[[
    'MAXENT(donor)',
    'MM(donor)',
    'WMM(donor)',
    'MAXENT(acceptor)',
    'MM(acceptor)',
    'WMM(acceptor)',
]].astype(float)


# conservation scores

conservation_scores_df = df[[
    'phyloP(acceptor_in)',
    'phyloP(acceptor_out)',
    'phyloP(donor_in)',
    'phyloP(donor_out)',
    'phastCons(acceptor_in)',
    'phastCons(acceptor_out)',
    'phastCons(donor_in)',
    'phastCons(donor_out)',
]].fillna('0').astype(float)


# circRNAs databases: num_db

num_db_df = df[[
    'num_db'
]]


# G score
G_score_df = df[[
    'G score(d=10)',
    'G score(d=5)'
]]


scoring_table_df = pd.concat(
    [
        base_df,
        num_reads_df,
        RPM_df,
        MCS_df,
        TSI_df,
        evidence_df,
        evidence_num_plus_df,
        flanking_regions_df,
        flanking_regions_df_2,
        splicing_scores_df,
        conservation_scores_df,
        num_db_df,
        G_score_df
    ],
    axis=1
)


scoring_table_df.reset_index().to_csv(snakemake.output[0], sep='\t', index=False)
