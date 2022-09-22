import pandas as pd
from operator import itemgetter


event_id_df = pd.read_csv(
    snakemake.input[0],
    sep='\t',
    dtype='object',
    usecols=[4]
).set_index('event_id')


dfs = [
    pd.read_csv(file_, sep='\t', dtype='object')
    for file_ in snakemake.input[1:]
]

totalRNA_read_count_dfs = [
    df[['event_id', 'totalRNA']].set_index('event_id')
    for df in dfs
]

merged_df = pd.concat([event_id_df] + totalRNA_read_count_dfs, axis=1)


sample_ids = [
    f'{tissue}({dataset_id})'
    for dataset_id, species, tissue in snakemake.params.groups
]

sample_id_mapper = snakemake.params.sample_id_mapper

new_sample_infos = [sample_id_mapper[id_] for id_ in sample_ids]
new_sample_ids, new_sample_order = zip(*new_sample_infos)

merged_df.columns = new_sample_ids


sorted_new_sample_ids, _ = zip(*sorted(new_sample_infos, key=itemgetter(1)))
sorted_new_sample_ids = list(sorted_new_sample_ids)

merged_df = merged_df[sorted_new_sample_ids]


merged_df.reset_index().to_csv(snakemake.output[0], sep='\t', index=False)