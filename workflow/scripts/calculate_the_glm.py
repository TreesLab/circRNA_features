import pandas as pd
import statsmodels.api as sm


df = pd.read_csv(snakemake.input[0], sep='\t').set_index('event_id')

data_endog = df.iloc[:, 0]
data_exog = df.iloc[:, 1:]


# GLM (full model)
binom_model = sm.GLM(data_endog, sm.add_constant(data_exog), family=sm.families.Binomial())
binom_results = binom_model.fit(use_t=True)


# GLM (reduced model)
reduced_data_exog = [data_exog.drop(col, axis=1) for col in data_exog.columns]

reduced_binom_results = []
for _data in reduced_data_exog:
    reduced_binom_model = sm.GLM(data_endog, sm.add_constant(_data), family=sm.families.Binomial())
    reduced_binom_results.append(reduced_binom_model.fit(use_t=True))


# GLM (specific reduced model)
specific_reduced_data_exog = []
specific_reduced_data_exog.append(
    data_exog.drop(
        [
            '#BSJ_reads(totalRNA)',
            'Type length(mature)',
            'MCS-detail(sample)'
        ],
        axis=1
    )
)

specific_reduced_binom_results = []
for _data in specific_reduced_data_exog:
    specific_reduced_binom_model = sm.GLM(data_endog, sm.add_constant(_data), family=sm.families.Binomial())
    specific_reduced_binom_results.append(specific_reduced_binom_model.fit(use_t=True))



# output results
out_data = []
out_data += [binom_results.pseudo_rsquared()]
out_data += binom_results.params.values.tolist()
out_data += binom_results.tvalues.values.tolist()
out_data += binom_results.bse.values.tolist()
out_data += binom_results.pvalues.values.tolist()
out_data += [result.pseudo_rsquared() for result in reduced_binom_results]
out_data += [result.pseudo_rsquared() for result in specific_reduced_binom_results]

out_data_titles = []
out_data_titles += ['r2_all']
out_data_titles += [f'beta{i}' for i in range(1 + data_exog.shape[1])]
out_data_titles += [f't{i}' for i in range(1 + data_exog.shape[1])]
out_data_titles += [f'sd{i}' for i in range(1 + data_exog.shape[1])]
out_data_titles += [f'p{i}' for i in range(1 + data_exog.shape[1])]
out_data_titles += [f'r2_all-feature{i}' for i in range(1, 1 + data_exog.shape[1])]
out_data_titles += ['r2_all-f1-f3-f4']

our_results = pd.Series(
    out_data,
    index=out_data_titles,
    name=f'{snakemake.wildcards.tissue}({snakemake.wildcards.dataset_id})'
)

our_results.to_csv(snakemake.output[0], sep='\t')
