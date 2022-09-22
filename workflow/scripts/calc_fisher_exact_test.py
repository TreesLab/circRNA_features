import pandas as pd
from scipy.stats import fisher_exact


df = pd.read_csv(snakemake.input[0], sep='\t', index_col=0, header=[0, 1])

a, b, c, d = (f'{name1}({name2})' for name1, name2 in df.columns)
frac_index = [
    (
        'a/(a+b)',
        f'a: {a}, b: {b}'
    ),
    (
        'c/(c+d)',
        f'c: {c}, d: {d}'
    )
]

frac_df = df.apply(
    lambda s: pd.Series([float(s[0]) / (s[0] + s[1]), float(s[2]) / (s[2] + s[3])], index=frac_index),
    axis=1
)

results_df = df.apply(
    lambda s: pd.Series(fisher_exact(s.values.reshape((2, 2))), index=[('odds_ratio', ' '), ('p_value', 'two-sided')]),
    axis=1
)

results_df_less = df.apply(
    lambda s: pd.Series(fisher_exact(s.values.reshape((2, 2)), alternative='less')[1], index=[('p_value', 'less')]),
    axis=1
)

results_df_greater = df.apply(
    lambda s: pd.Series(fisher_exact(s.values.reshape((2, 2)), alternative='greater')[1], index=[('p_value', 'greater')]),
    axis=1
)

df_with_results = pd.concat([df, frac_df, results_df, results_df_less, results_df_greater], axis=1)
df_with_results.to_csv(snakemake.output[0], sep='\t')
