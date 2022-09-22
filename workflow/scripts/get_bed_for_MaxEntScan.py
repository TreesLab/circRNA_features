import pandas as pd


circ_df = pd.read_csv(snakemake.input[0], sep='\t', dtype='object').astype({'donor': int, 'acceptor': int})


circ_df_p = circ_df[circ_df['strand'] == '+']
circ_df_n = circ_df[circ_df['strand'] == '-']


circ_df_p_d = circ_df_p.assign(
    pos1=circ_df_p['donor'] - 3,
    pos2=circ_df_p['donor'] + 6,
    score='0'
)[['chr', 'pos1', 'pos2', 'event_id', 'score', 'strand']]

circ_df_n_d = circ_df_n.assign(
    pos1=circ_df_n['donor'] - 7,
    pos2=circ_df_n['donor'] + 2,
    score='0'
)[['chr', 'pos1', 'pos2', 'event_id', 'score', 'strand']]

circ_df_n_a = circ_df_n.assign(
    pos1=circ_df_n['acceptor'] - 3,
    pos2=circ_df_n['acceptor'] + 20,
    score='0'
)[['chr', 'pos1', 'pos2', 'event_id', 'score', 'strand']]

circ_df_p_a = circ_df_p.assign(
    pos1=circ_df_p['acceptor'] - 21,
    pos2=circ_df_p['acceptor'] + 2,
    score='0'
)[['chr', 'pos1', 'pos2', 'event_id', 'score', 'strand']]


circ_df_d = pd.concat([circ_df_p_d, circ_df_n_d]).sort_index()
circ_df_a = pd.concat([circ_df_p_a, circ_df_n_a]).sort_index()

circ_df_d.to_csv(snakemake.output.five_prime_bed, sep='\t', index=False, header=False)
circ_df_a.to_csv(snakemake.output.three_prime_bed, sep='\t', index=False, header=False)
