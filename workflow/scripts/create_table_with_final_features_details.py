import pandas as pd
import re
from operator import itemgetter



circ_df = pd.read_csv(snakemake.input.circRNAs, sep='\t', dtype='object')
event_id_df = circ_df[['event_id']]


def append_all_features(circ_df, *feature_dfs):
    for feature_df in feature_dfs:
        circ_df = circ_df.merge(feature_df, on='event_id', how='left')

    return circ_df


# check_ambiguous
check_ambiguous_df = pd.read_csv(
    snakemake.input.check_ambiguous,
    sep='\t'
).rename(
    {
        'circRNA': 'event_id'
    },
    axis=1
).set_index(
    'event_id'
).assign(
    has_ambiguity=lambda df: (df.sum(axis=1) > 0).apply(int)
)[[
    'has_ambiguity'
]]


# circAtlas table
circAtlas_df = pd.read_csv(
    snakemake.input.circAtlas,
    sep='\t',
    dtype='object'
)

# circAtlas table: algorithms
algorithms_df = circAtlas_df[[
    'event_id',
    'Algorithm(CIRI2)',
    'Algorithm(DCC)',
    'Algorithm(find_circ)',
    'Algorithm(CIRCexplorer2)'
]].set_index(
    'event_id'
).astype(int).assign(
    algorithms=lambda df: df.apply(sum, axis=1)
).assign(
    **{
        'algorithms > 1': lambda df: (df['algorithms'] > 1).apply(int),
    }
)[[
    'algorithms > 1'
]]

# circAtlas table: Type length(mature)
type_length_mature_df = circAtlas_df[[
    'event_id',
    'Type length(mature)'
]].set_index('event_id')


# circAtlas table: MCS
MCS_df = circAtlas_df[[
    'event_id',
    'MCS-detail(species)',
    'MCS-detail(tissue)',
    'MCS-detail(sample)',
    'Multiple Conservation Score(MCS)',
    'Tissue Specifity Index'
]].set_index('event_id')



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
).drop_duplicates().merge(
    event_id_df,
    on='event_id',
    how='outer'
).fillna('0')



# conservation scores
conservation_scores_df = pd.read_csv(
    snakemake.input.conservation_scores,
    sep='\t'
).set_index('event_id')



# donor & acceptor: annotation check
check_annotated_df = pd.read_csv(
    snakemake.input.check_annotated,
    sep='\t',
    dtype='object',
    usecols=[4, 6, 7, 8]
).rename(
    {
        'circ_id': 'event_id'
    },
    axis=1
).assign(
    both=lambda df: ((df.iloc[:, 1] == '1') & (df.iloc[:, 2] == '1')).apply(int).apply(str)
).rename(
    {
        'both': 'both_donor_acceptor_are_at_annotated_boundary'
    },
    axis=1
).set_index('event_id')


# donor & acceptor: AS events check
check_AS_df = pd.read_csv(
    snakemake.input.check_AS_events,
    sep='\t',
    dtype='object',
    usecols=[4, 6, 8]
).assign(
    both=lambda df: ((df.iloc[:, 1] == '1') & (df.iloc[:, 2] == '1')).apply(int).apply(str)
).rename(
    {
        'both': 'both_donor_acceptor_have_AS_events'
    },
    axis=1
).set_index('event_id')



# RCS
RCS_df = pd.read_csv(
    snakemake.input.RCS_pairs,
    sep='\t',
    dtype='object',
    usecols=[0, 1, 5, 6]
).rename(
    {
        'circRNA_id': 'event_id'
    },
    axis=1
).merge(event_id_df, on='event_id', how='right').set_index(
    'event_id'
).fillna('0').astype(int).assign(
    **{
        '#RCS_across > 0': lambda df: (df['#RCS_across_flanking_introns'] > 0).apply(int),
        '#RCS_across - #RCS_within > 0': lambda df: ((df['#RCS_across_flanking_introns'] - df['#RCS_within_donor\'s_intron'] - df['#RCS_within_acceptor\'s_intron']) > 0).apply(int)
    }
).drop(
    [
        '#RCS_across_flanking_introns',
        '#RCS_within_donor\'s_intron',
        '#RCS_within_acceptor\'s_intron'
    ],
    axis=1
)



# flanking regions: RBP pairs on flanking 1kb
RBP_df = pd.read_csv(
    snakemake.input.RBP_pairs,
    sep='\t',
    dtype='object',
    names=['event_id', '#RBP_pairs_on_flanking_1k'],
    usecols=[0, 1]
).set_index('event_id').assign(
    has_common_RBPs_on_flanking_1k=lambda df: (df['#RBP_pairs_on_flanking_1k'].apply(int) > 0).apply(int)
)[[
    'has_common_RBPs_on_flanking_1k'
]]



# donor & acceptor: splicing scores
splicing_scores_df = pd.read_csv(
    snakemake.input.splicing_scores,
    sep='\t',
    dtype='object'
)[[
    'event_id',
    'MAXENT(acceptor)',
    'MAXENT(donor)',
    'MM(acceptor)',
    'MM(donor)',
    'WMM(acceptor)',
    'WMM(donor)',
]].set_index('event_id')



# transCirc table: evidences
evidences_df = pd.read_csv(
    snakemake.input.transCirc,
    sep='\t',
    dtype='object'
)[[
    'event_id',
    'evidences_num',
    'evidences_score'
]].merge(event_id_df, on='event_id', how='right').set_index(
    'event_id'
).fillna('0')


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
).merge(event_id_df, on='event_id', how='right').set_index(
    'event_id'
)[[
    'has_miRNA_binding_sites_across_junction'
]].fillna(0)


# circRNA junction: RBP-binding sites across the circRNA junction
cross_junc_RBPs_high_cons_df = pd.read_csv(
    snakemake.input.cross_junc_RBPs_high_cons,
    sep='\t',
    dtype='object',
    usecols=[0, 1]
).rename(
    {
        'seq_id': 'event_id',
        'count': '#RBP_with_binding_sites_across_junction_high_cons'
    },
    axis=1
).set_index('event_id').assign(
    has_RBP_binding_sites_across_junction_high_cons=lambda df: (df['#RBP_with_binding_sites_across_junction_high_cons'].apply(int) > 0).apply(int)
).merge(event_id_df, on='event_id', how='right').set_index(
    'event_id'
)[[
    'has_RBP_binding_sites_across_junction_high_cons'
]].fillna(0)


# evidence_num + "miRNA/RBP cross junction" 
evidence_num_plus_df = pd.concat(
    [
        evidences_df[['evidences_num']].astype(int),
        cross_junc_miRNAs_df,
        cross_junc_RBPs_high_cons_df
    ],
    axis=1
).assign(
    evidence_num_plus=lambda df: df.sum(axis=1).astype(int)
)[[
    'evidence_num_plus'
]]

# curated_circRNAs_db
curated_circRNAs_db = pd.read_csv(
    snakemake.input.curated_circRNAs_db,
    sep='\t',
    dtype='object'
)[[
    'event_id',
    'CircR2diseaseV2(well_confirmed)',
    'RT-independent'
]]


# circRNAs databases
db_df = pd.read_csv(
    snakemake.input.circRNAs_db,
    sep='\t',
    dtype='object'
)


# sample features
sample_dfs = [
    pd.read_csv(file_, sep='\t', dtype='object', usecols=[0, 1, 2, 3]).set_index('event_id')
    for file_ in snakemake.input.sample_features
]

sample_ids = [
    f'{tissue}({dataset_id})'
    for dataset_id, species, tissue in snakemake.params.groups
]

sample_id_mapper = snakemake.params.sample_id_mapper

new_sample_infos = [sample_id_mapper[id_] for id_ in sample_ids]
new_sample_ids, new_sample_order = zip(*new_sample_infos)
new_sample_ids = [re.sub(r'\(.+\)$', '', id_) for id_ in new_sample_ids]


sample_data = sorted(zip(sample_dfs, new_sample_ids, new_sample_order), key=itemgetter(2))

for df, id_, _ in sample_data:
    df.columns = [f'{col}({id_})' for col in df.columns]

sorted_sample_dfs, _, _ = list(zip(*sample_data))
sorted_sample_dfs = list(sorted_sample_dfs)

sample_features_df = append_all_features(
    event_id_df,
    *sorted_sample_dfs
)


# merge all features
merged_df = append_all_features(
    event_id_df,
    check_ambiguous_df,
    algorithms_df,
    type_length_mature_df,
    circFLseq_df,
    MCS_df,
    conservation_scores_df,
    check_annotated_df,
    check_AS_df,
    RCS_df,
    RBP_df,
    splicing_scores_df,
    evidence_num_plus_df,
    evidences_df,
    curated_circRNAs_db,
    db_df
)

merged_df = merged_df[merged_df['has_ambiguity'] == 0].drop('has_ambiguity', axis=1)

merged_df = merged_df[[
    'event_id',
    'algorithms > 1',
    'Type length(mature)',
    'is_in_circFLseq',
    'MCS-detail(species)',
    'MCS-detail(tissue)',
    'MCS-detail(sample)',
    'Multiple Conservation Score(MCS)',
    'Tissue Specifity Index',
    'phyloP(acceptor_in)',
    'phyloP(acceptor_out)',
    'phyloP(donor_in)',
    'phyloP(donor_out)',
    'phastCons(acceptor_in)',
    'phastCons(acceptor_out)',
    'phastCons(donor_in)',
    'phastCons(donor_out)',
    'donor_site_at_the_annotated_boundary',
    'acceptor_site_at_the_annotated_boundary',
    'both_donor_acceptor_are_at_annotated_boundary',
    'donor_acceptor_sites_at_the_same_transcript_isoform',
    'has_AS_event(donor)',
    'has_AS_event(acceptor)',
    'both_donor_acceptor_have_AS_events',
    'MAXENT(acceptor)',
    'MAXENT(donor)',
    'MM(acceptor)',
    'MM(donor)',
    'WMM(acceptor)',
    'WMM(donor)',
    '#RCS_across > 0',
    '#RCS_across - #RCS_within > 0',
    'has_common_RBPs_on_flanking_1k',
    'evidence_num_plus',
    'evidences_score',
    'CircR2diseaseV2(well_confirmed)',
    'RT-independent',
    'CIRCpedia_v2',
    'CSCD_v2',
    'CircRic',
    'MiOncoCirc',
    'TSCD',
    'circBase',
    'circRNADb',
    'exoRBase',
    'num_db',
]].rename(
    {
        'algorithms > 1': 'Detected by multiple tools (algorithms >=2) (Yes: 1; No: 0)',
        'Type length(mature)': 'Full-length circular sequence (by CircAtlas-CIRI-full) (Yes: 1; No: 0)',
        'is_in_circFLseq': 'Full-length circular sequence (by circFLseq) (Yes: 1; No: 0)',
        'MCS-detail(species)': '#species (CircAtlas MCS-detail)',
        'MCS-detail(tissue)': '#tissue (CircAtlas MCS-detail)',
        'MCS-detail(sample)': '#sample (CircAtlas MCS-detail)',
        'Multiple Conservation Score(MCS)': 'Multiple Conservation Score(MCS; CircAtlas )',
        'Tissue Specifity Index': 'Tissue Specifity Index (CircAtlas)',
        'phyloP(acceptor_in)': 'phyloP (acceptor_in exon)',
        'phyloP(acceptor_out)': 'phyloP (acceptor_in intron)',
        'phyloP(donor_in)': 'phyloP (donor_in exon)',
        'phyloP(donor_out)': 'phyloP (donor_in intron)',
        'phastCons(acceptor_in)': 'phastCons (acceptor_in exon)',
        'phastCons(acceptor_out)': 'phastCons (acceptor_in intron)',
        'phastCons(donor_in)': 'phastCons (donor_in exon)',
        'phastCons(donor_out)': 'phastCons (donor_in intron)',
        'acceptor_site_at_the_annotated_boundary': 'acceptor_site_at_the_annotated_boundary (Yes: 1; No: 0)',
        'donor_site_at_the_annotated_boundary': 'donor_site_at_the_annotated_boundary (Yes: 1; No: 0)',
        'both_donor_acceptor_are_at_annotated_boundary': 'both_donor_acceptor_at_annotated_boundary (Yes: 1; No: 0)',
        'donor_acceptor_sites_at_the_same_transcript_isoform': 'both_donor_acceptor_sites_at_the_same_transcript_isoform (Yes: 1; No: 0)',
        'has_AS_event(acceptor)': 'Acceptor site undergoing AS (Yes: 1; No: 0)',
        'has_AS_event(donor)': 'donor site undergoing AS (Yes: 1; No: 0)',
        'both_donor_acceptor_have_AS_events': 'both_donor_acceptor_undergoing AS (Yes: 1; No: 0)',
        '#RCS_across > 0': '#RCS_across > 0 (Yes: 1; No: 0)',
        '#RCS_across - #RCS_within > 0': '#RCS_across - #RCS_within > 0 (Yes: 1; No: 0)',
        'has_common_RBPs_on_flanking_1k': 'has_common_RBPs_on_flanking_1k (based on ENCORI) (Yes: 1; No: 0)',
        'MAXENT(acceptor)': 'MEM (acceptor)',
        'MAXENT(donor)': 'MEM (donor)',
        'MM(acceptor)': 'FMM (acceptor)',
        'MM(donor)': 'FMM (donor)',
        'WMM(acceptor)': 'WMM (acceptor)',
        'WMM(donor)': 'WMM (donor)',
        'evidence_num_plus': 'evidence_num_TransCirc+miRNA/RBP sites across the BSJ',
        'evidences_score': 'evidences_score (translational potential)',
    },
    axis=1
)


merged_df = append_all_features(
    merged_df,
    sample_features_df
)


# output
merged_df.to_csv(snakemake.output[0], sep='\t', index=False)


