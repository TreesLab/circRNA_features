import pandas as pd
import numpy as np
import math
import subprocess as sp
from scipy.stats import ranksums
from itertools import product

import logging

logging.basicConfig(filename=snakemake.log[0], encoding='utf-8', level=logging.INFO)


df = pd.read_csv(snakemake.input[0], sep='\t').set_index('event_id')


not_detected_col = list(filter(lambda col: col.startswith('Not detected'), df.columns))[0]
not_detected_df = df[df[not_detected_col] == 1].iloc[:, 2:]
detected_df = df[df[not_detected_col] == 0].iloc[:, 2:]

not_depleted_col = list(filter(lambda col: col.startswith('Not depleted'), df.columns))[0]
not_depleted_df = df[df[not_depleted_col] == 1].iloc[:, 2:]
depleted_df = df[df[not_depleted_col] == 0].iloc[:, 2:]


def get_test_pairs(df1, df2):
    for df1_col, df2_col in zip(df1.iteritems(), df2.iteritems()):
        assert df1_col[0] == df2_col[0]
        yield df1_col[0], (df1_col[1], df2_col[1])


logging.debug(next(get_test_pairs(not_detected_df, detected_df)))
logging.debug(next(get_test_pairs(not_depleted_df, depleted_df)))


def r_ranksum(df1, df2):
    data_1 = ','.join(df1.astype(str))
    data_2 = ','.join(df2.astype(str))
    data = f'{data_1}\n{data_2}\n'

    with sp.Popen(
        [
            'Rscript',
            snakemake.params.r_ranksum
        ],
        stdin=sp.PIPE,
        stdout=sp.PIPE,
        encoding='utf-8'
    ) as p:

        res, err_msg = p.communicate(data)

        pvalue, conf_int_low, conf_int_high, differnce_in_location = res.split()

        return pvalue, conf_int_low, conf_int_high, differnce_in_location


not_detected_results = [
    (
        col,
        pairs[0].median(),
        pairs[1].median(),
        ranksums(*pairs).pvalue,
        ranksums(*pairs, alternative='less').pvalue,
        ranksums(*pairs, alternative='greater').pvalue,
        ranksums(*pairs).statistic / math.sqrt(len(pairs[0]) + len(pairs[1])),
        np.median([v1 - v2 for v1, v2 in product(*pairs)]),
        *r_ranksum(*pairs)
    )
    for col, pairs in get_test_pairs(not_detected_df, detected_df)
]

not_detected_results_df = pd.DataFrame(
    not_detected_results,
    columns=[
        'features',
        'median(Not detected = 1)',
        'median(Not detected = 0)',
        'p_value(two-sided)(ranksums)',
        'p_value(less)(ranksums)',
        'p_value(greater)(ranksums)',
        'effect_size(two-sided)(ranksums)',
        'differnce_in_location(ranksums)',
        'p_value(two-sided)(R_ranksums)',
        'conf_int_low(two-sided)(R_ranksums)',
        'conf_int_high(two-sided)(R_ranksums)',
        'differnce_in_location(R_ranksums)',
    ]
)


not_depleted_results = [
    (
        col,
        pairs[0].median(),
        pairs[1].median(),
        ranksums(*pairs).pvalue,
        ranksums(*pairs, alternative='less').pvalue,
        ranksums(*pairs, alternative='greater').pvalue,
        ranksums(*pairs).statistic / math.sqrt(len(pairs[0]) + len(pairs[1])),
        np.median([v1 - v2 for v1, v2 in product(*pairs)]),
        *r_ranksum(*pairs)
    )
    for col, pairs in get_test_pairs(not_depleted_df, depleted_df)
]

not_depleted_results_df = pd.DataFrame(
    not_depleted_results,
    columns=[
        'features',
        'median(Not depleted = 1)',
        'median(Not depleted = 0)',
        'p_value(two-sided)(ranksums)',
        'p_value(less)(ranksums)',
        'p_value(greater)(ranksums)',
        'effect_size(two-sided)(ranksums)',
        'differnce_in_location(ranksums)',
        'p_value(two-sided)(R_ranksums)',
        'conf_int_low(two-sided)(R_ranksums)',
        'conf_int_high(two-sided)(R_ranksums)',
        'differnce_in_location(R_ranksums)',
    ]
)

not_detected_results_df.to_csv(snakemake.output[0], sep='\t', index=False)
not_depleted_results_df.to_csv(snakemake.output[1], sep='\t', index=False)
