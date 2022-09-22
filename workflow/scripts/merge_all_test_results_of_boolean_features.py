import pandas as pd
from operator import itemgetter


import logging

logging.basicConfig(filename=snakemake.log[0], encoding='utf-8', level=logging.DEBUG)


def get_sample_info_from_file_path(file_path):
    split_path = file_path.split('/')
    logging.debug(split_path)

    dataset_id = split_path[2]
    species = split_path[3]
    tissue = split_path[4]

    return dataset_id, species, tissue


all_samples_info = [
    get_sample_info_from_file_path(file_path)
    for file_path in snakemake.input
]

all_sample_headers = [
    f"{tissue}({dataset_id})"
    for dataset_id, species, tissue in all_samples_info
]


def create_new_column_names(old_column_names, sample_name):
    new_column_names = ['{}({})'.format(*names) if names[1] != ' ' else names[0] for names in old_column_names]
    new_column_names = pd.MultiIndex.from_product([[sample_name], new_column_names])
    return new_column_names


all_sample_dfs = []
for file_path, sample_header in zip(snakemake.input, all_sample_headers):
    df = pd.read_csv(file_path, sep='\t', dtype='object', index_col=0, header=[0, 1])
    df.columns = create_new_column_names(df.columns, sample_header)
    all_sample_dfs.append(df.stack())

merged_df = pd.concat(
    all_sample_dfs,
    axis=1
)


# modify headers

sample_id_mapper = snakemake.params.sample_id_mapper

column_names = merged_df.columns

new_sample_infos = [sample_id_mapper[col] for col in column_names]
new_sample_ids, new_sample_order = zip(*new_sample_infos)

merged_df.columns = new_sample_ids

sorted_new_sample_ids, _ = zip(*sorted(new_sample_infos, key=itemgetter(1)))
sorted_new_sample_ids = list(sorted_new_sample_ids)

merged_df = merged_df[sorted_new_sample_ids]




merged_df.to_csv(snakemake.output[0], sep='\t')



