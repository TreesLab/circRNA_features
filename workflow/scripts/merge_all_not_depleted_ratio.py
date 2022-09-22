import pandas as pd
from operator import itemgetter


df = pd.DataFrame([], columns=['event_id'])

for i, file_ in enumerate(snakemake.input.files, start=1):
    ratio_df = pd.read_csv(
        file_,
        sep='\t',
        dtype='object'
    )[[
        'event_id',
        'not_depleted_ratio'
    ]].rename(
        {
            'not_depleted_ratio': f'not_depleted_ratio_{i}'
        },
        axis=1
    )

    df = df.merge(ratio_df, on='event_id', how='outer')


sample_ids = [
    f'{tissue}({dataset_id})'
    for dataset_id, species, tissue in snakemake.params.groups
]

sample_id_mapper = snakemake.params.sample_id_mapper

new_sample_infos = [sample_id_mapper[id_] for id_ in sample_ids]
new_sample_ids, new_sample_order = zip(*new_sample_infos)


column_names = [
    f'not_depleted_ratio ({id_})'
    for id_ in new_sample_ids
]

df.columns = ['event_id'] + column_names


sorted_columns, _ = zip(*sorted(zip(column_names, new_sample_order), key=itemgetter(1)))

df = df[['event_id'] + list(sorted_columns)]


df.to_csv(snakemake.output[0], sep='\t', index=False)
