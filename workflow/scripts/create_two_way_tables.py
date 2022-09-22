import pandas as pd
import logging


logging.basicConfig(filename=snakemake.log[0], encoding='utf-8', level=logging.DEBUG)

logging.info('loading features table')
df = pd.read_csv(snakemake.input[0], sep='\t', dtype='object', index_col=0)

logging.info('generating "not detected" and "not depleted" columns')
column_names = df.columns

not_detected_col = column_names[0]
not_depleted_col = column_names[1]
bool_features = column_names[2:]

not_detected_table_columns = pd.MultiIndex.from_product(
    [
        ('Not detected', 'detected'),
        ('yes', 'no')
    ],
    names=['', 'has_feature_or_not']
)

not_depleted_table_columns = pd.MultiIndex.from_product(
    [
        ('Not depleted', 'depleted'),
        ('yes', 'no')
    ],
    names=['', 'has_feature_or_not']
)

table_index = pd.Index(bool_features, name='features')


def create_two_way_table(df, feature1, feature2):
    two_way_df = df[
        [
            feature1,
            feature2
        ]
    ].assign(
        count=1
    ).groupby(
        [
            feature1,
            feature2
        ]
    ).agg('count')

    logging.debug(two_way_df)

    two_way_dict = two_way_df.to_dict()['count']

    logging.debug(two_way_dict)

    results = [
        two_way_dict.get(('1', '1'), 0),
        two_way_dict.get(('1', '0'), 0),
        two_way_dict.get(('0', '1'), 0),
        two_way_dict.get(('0', '0'), 0)
    ]

    return results

logging.info('creating two way tables with "not detected" column')
not_detected_table = pd.DataFrame(
    (
        create_two_way_table(df, not_detected_col, feature)
        for feature in bool_features
    ),
    index=table_index,
    columns=not_detected_table_columns
)

not_detected_table.to_csv(snakemake.output.not_detected_table, sep='\t')


logging.info('creating two way tables with "not depleted" column')
not_depleted_table = pd.DataFrame(
    (
        create_two_way_table(df, not_depleted_col, feature)
        for feature in bool_features
    ),
    index=table_index,
    columns=not_depleted_table_columns
)

not_depleted_table.to_csv(snakemake.output.not_depleted_table, sep='\t')


logging.info('all done!')
