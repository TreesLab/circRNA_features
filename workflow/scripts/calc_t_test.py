import pandas as pd
from scipy.stats import ttest_ind

import logging

logging.basicConfig(filename=snakemake.log[0], encoding='utf-8', level=logging.DEBUG)


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


not_detected_results = [
    (
        col,
        pairs[0].mean(),
        pairs[1].mean(),
        ttest_ind(*pairs, nan_policy='raise').pvalue,
        ttest_ind(*pairs, nan_policy='raise', alternative='less').pvalue,
        ttest_ind(*pairs, nan_policy='raise', alternative='greater').pvalue
    )
    for col, pairs in get_test_pairs(not_detected_df, detected_df)
]

not_detected_results_df = pd.DataFrame(
    not_detected_results,
    columns=[
        'features',
        'mean(Not detected = 1)',
        'mean(Not detected = 0)',
        'p_value(two-sided)(t_test)',
        'p_value(less)(t_test)',
        'p_value(greater)(t_test)'
    ]
)


not_depleted_results = [
    (
        col,
        pairs[0].mean(),
        pairs[1].mean(),
        ttest_ind(*pairs, nan_policy='raise').pvalue,
        ttest_ind(*pairs, nan_policy='raise', alternative='less').pvalue,
        ttest_ind(*pairs, nan_policy='raise', alternative='greater').pvalue
    )
    for col, pairs in get_test_pairs(not_depleted_df, depleted_df)
]

not_depleted_results_df = pd.DataFrame(
    not_depleted_results,
    columns=[
        'features',
        'mean(Not depleted = 1)',
        'mean(Not depleted = 0)',
        'p_value(two-sided)(t_test)',
        'p_value(less)(t_test)',
        'p_value(greater)(t_test)'
    ]
)

not_detected_results_df.to_csv(snakemake.output[0], sep='\t', index=False)
not_depleted_results_df.to_csv(snakemake.output[1], sep='\t', index=False)
