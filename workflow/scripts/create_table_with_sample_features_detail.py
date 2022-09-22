import pandas as pd


rm_amb_df = pd.read_csv(snakemake.input.rm_amb, sep='\t', dtype='object')
boolean_df = pd.read_csv(snakemake.input.boolean, sep='\t', dtype='object')
scoring_df = pd.read_csv(snakemake.input.scoring, sep='\t', dtype='object')


treatment_columns = snakemake.params.treatments
RPM_columns = list(filter(lambda col: col.startswith('RPM_'), rm_amb_df.columns))

totalRNA_RPM_col = list(filter(lambda col: col.startswith(f'RPM_totalRNA'), RPM_columns))[0]

RNaseR_polyAminus_col = list(filter(lambda col: col in ['RNaseR', 'polyA_minus'], treatment_columns))[0]
RNaseR_polyAminus_RPM_col = list(filter(lambda col: col.startswith(f'RPM_{RNaseR_polyAminus_col}'), RPM_columns))[0]


base_df = rm_amb_df[[
    'event_id',
    'totalRNA',
    RNaseR_polyAminus_col,
    totalRNA_RPM_col,
    RNaseR_polyAminus_RPM_col
]]


not_depleted_df = pd.read_csv(
    snakemake.input.not_depleted,
    sep='\t',
    dtype='object'
)[[
    'not_depleted_ratio'
]]


wanted_boolean_df = boolean_df[[
    'algorithms > 1',
    'Type length(mature)',
    'is_in_circFLseq',
    'acceptor_site_at_the_annotated_boundary',
    'donor_site_at_the_annotated_boundary',
    'both_donor_acceptor_are_at_annotated_boundary',
    'donor_acceptor_sites_at_the_same_transcript_isoform',
    'has_AS_event(acceptor)',
    'has_AS_event(donor)',
    'both_donor_acceptor_have_AS_events',
    '#RCS_across > 0',
    '#RCS_across - #RCS_within > 0',
    'has_common_RBPs_on_flanking_1k',

]]


wanted_scoring_df = scoring_df[[
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
    'MAXENT(acceptor)',
    'MAXENT(donor)',
    'MM(acceptor)',
    'MM(donor)',
    'WMM(acceptor)',
    'WMM(donor)',
    'evidence_num_plus',
    'evidences_num',
    'evidences_score',
]]



merged_df = pd.concat(
    [
        base_df,
        not_depleted_df,
        wanted_boolean_df,
        wanted_scoring_df
    ],
    axis=1
)


merged_df = merged_df[[
    'event_id',
    'totalRNA',
    RNaseR_polyAminus_col,
    'not_depleted_ratio',
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
    'acceptor_site_at_the_annotated_boundary',
    'donor_site_at_the_annotated_boundary',
    'both_donor_acceptor_are_at_annotated_boundary',
    'donor_acceptor_sites_at_the_same_transcript_isoform',
    'has_AS_event(acceptor)',
    'has_AS_event(donor)',
    'both_donor_acceptor_have_AS_events',
    '#RCS_across > 0',
    '#RCS_across - #RCS_within > 0',
    'has_common_RBPs_on_flanking_1k',
    'MAXENT(acceptor)',
    'MAXENT(donor)',
    'MM(acceptor)',
    'MM(donor)',
    'WMM(acceptor)',
    'WMM(donor)',
    'evidence_num_plus',
    'evidences_num',
    'evidences_score',
]].rename(
    {
        'totalRNA': '#BSJ_reads(mock)',
        RNaseR_polyAminus_col: '#BSJ_reads(treated)',
        'not_depleted_ratio': 'treated/mock ratio',
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
        'donor_acceptor_sites_at_the_same_transcript_isoform': 'donor_acceptor_sites_at_the_same_transcript_isoform (Yes: 1; No: 0)',
        'has_AS_event(acceptor)': 'Acceptor site undergoing AS (Yes: 1; No: 0)',
        'has_AS_event(donor)': 'donor site undergoing AS (Yes: 1; No: 0)',
        'both_donor_acceptor_have_AS_events': 'both_donor_acceptor_undergoing AS (Yes: 1; No: 0)',
        '#RCS_across > 0': '#RCS_across > 0 (Yes: 1; No: 0)',
        '#RCS_across - #RCS_within > 0': '#RCS_across - #RCS_within > 0 (Yes: 1; No: 0)',
        'has_common_RBPs_on_flanking_1k': 'has_common_RBPs_on_flanking_1k (based on ENCORI) (Yes: 1; No: 0)',
        'MAXENT(acceptor)': 'MAXENT (acceptor)',
        'MAXENT(donor)': 'MAXENT (donor)',
        'MM(acceptor)': 'MM (acceptor)',
        'MM(donor)': 'MM (donor)',
        'WMM(acceptor)': 'WMM (acceptor)',
        'WMM(donor)': 'WMM (donor)',
        'evidence_num_plus': 'evidence_num_TransCirc+miRNA/RBP sites across the BSJ',
        'evidences_num': 'evidences_num (translational potential)',
        'evidences_score': 'evidences_score (translational potential)',
    },
    axis=1
)


merged_df.to_csv(snakemake.output[0], sep='\t', index=False)

