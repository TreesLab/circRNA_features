import pandas as pd



def get_column_name(columns, text_starts):
    col_name = list(filter(lambda col: col.startswith(text_starts), columns))

    if col_name:
        return col_name[0]
    else:
        return None


def get_summary(df):
    columns = df.columns

    not_detected_column = get_column_name(columns, 'Not detected')
    not_depleted_column = get_column_name(columns, 'Not depleted')

    needed_columns = [
        'count',
        not_detected_column,
        not_depleted_column
    ]

    count_result = df.assign(count=1)[needed_columns].astype(int).sum().values.tolist()

    return count_result


summary_df = pd.DataFrame(
    (
        [file_.lstrip(snakemake.params.basename)] + get_summary(pd.read_csv(file_, sep='\t', dtype='object'))
        for file_ in snakemake.input
    ),
    columns=['filename', 'count', 'Not detected', 'Not depleted']
).assign(
    detected_and_depleted=lambda df: df['count'] - df['Not detected'] - df['Not depleted']
)[[
    'filename',
    'count',
    'Not detected',
    'Not depleted',
    'detected_and_depleted',
]].rename(
    {
        'count': '#circRNAs',
        'Not detected': '#(Not detected = 1)',
        'Not depleted': '#(Not depleted = 1)',
        'detected_and_depleted': '#(Not detected = 0 & Not depleted = 0)',
    },
    axis=1
).sort_values('filename').assign(
    ratio=lambda df: df['#(Not depleted = 1)'].astype(float) / df['#circRNAs']
).rename(
    {
        'ratio': '#(Not depleted = 1) / #circRNAs'    
    },
    axis=1
)

summary_df.to_csv(snakemake.output[0], sep='\t', index=False)
