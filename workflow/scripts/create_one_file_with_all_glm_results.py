import pandas as pd
from operator import itemgetter


dfs = [pd.read_csv(file_, sep='\t', index_col=0, dtype='object') for file_ in snakemake.input]
merged_df = pd.concat(dfs, axis=1)

sample_id_mapper = snakemake.params.sample_id_mapper

column_names = merged_df.columns

new_sample_infos = [sample_id_mapper[col] for col in column_names]
new_sample_ids, new_sample_order = zip(*new_sample_infos)

merged_df.columns = new_sample_ids

sorted_new_sample_ids, _ = zip(*sorted(new_sample_infos, key=itemgetter(1)))
sorted_new_sample_ids = list(sorted_new_sample_ids)

merged_df = merged_df[sorted_new_sample_ids]

merged_df.to_csv(snakemake.output[0], sep='\t')