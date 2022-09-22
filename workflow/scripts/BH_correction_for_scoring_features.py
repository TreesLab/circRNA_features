import pandas as pd
from statsmodels.stats import multitest

import logging
logging.basicConfig(filename=snakemake.log[0], encoding='utf-8', level=logging.DEBUG)


df = pd.read_csv(snakemake.input[0], sep='\t', dtype='object', index_col=[0, 1])

pv_cols = [
    'p_value(greater)(ranksums)',
    'p_value(less)(ranksums)',
    'p_value(two-sided)(ranksums)',
    'p_value(greater)(t_test)',
    'p_value(less)(t_test)',
    'p_value(two-sided)(t_test)'
]

pv_df = df.loc[(slice(None), pv_cols), :]
bh_pv_df = pv_df.transform(
    lambda x: multitest.multipletests(x.astype(float), method='fdr_bh')[1],
    axis=1
).rename(
    lambda x: x + "(BH)",
    level=1
)

merged_df = pd.concat([df, bh_pv_df]).sort_index()
merged_df.to_csv(snakemake.output[0], sep='\t')
