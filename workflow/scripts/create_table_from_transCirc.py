import pandas as pd


transCirc_df = pd.read_csv(snakemake.input[0], sep='\t', dtype='object')#, usecols=[0, 1, 2, 3])
circAtlas_df = pd.read_csv(snakemake.input[1], sep='\t', dtype='object', usecols=[0, 1])

transCirc_df_with_event_id = circAtlas_df.merge(transCirc_df, on='circAtlas_id', how='left').drop('circAtlas_id', axis=1)

transCirc_df_with_event_id.to_csv(snakemake.output[0], sep='\t', index=False)
