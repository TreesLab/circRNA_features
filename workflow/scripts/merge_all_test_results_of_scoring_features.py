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


all_samples_info_ranksum = [
    get_sample_info_from_file_path(file_path)
    for file_path in snakemake.input.ranksum_results
]

all_samples_info_t_test = [
    get_sample_info_from_file_path(file_path)
    for file_path in snakemake.input.t_test_results
]

assert all_samples_info_ranksum == all_samples_info_t_test

all_samples_info = all_samples_info_ranksum


all_sample_headers_ranksum = [
    f"{tissue}({dataset_id})"
    for dataset_id, species, tissue in all_samples_info_ranksum
]

all_sample_headers_t_test = [
    f"{tissue}({dataset_id})"
    for dataset_id, species, tissue in all_samples_info_t_test
]

assert all_sample_headers_ranksum == all_sample_headers_t_test

all_sample_headers = all_sample_headers_ranksum


def create_new_column_names(old_column_names, sample_name):
    new_column_names = pd.MultiIndex.from_product([[sample_name], old_column_names])
    return new_column_names


all_sample_dfs = []
for file_path_1, file_path_2, sample_header in zip(snakemake.input.ranksum_results, snakemake.input.t_test_results, all_sample_headers):
    df_1 = pd.read_csv(file_path_1, sep='\t', dtype='object', index_col=0)
    df_2 = pd.read_csv(file_path_2, sep='\t', dtype='object', index_col=0)
    df = pd.concat([df_1, df_2], axis=1)
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
