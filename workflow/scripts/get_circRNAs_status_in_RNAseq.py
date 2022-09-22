import pandas as pd


df = pd.read_csv(snakemake.input[0], sep='\t', dtype='object')

treatment_columns = snakemake.params.treatments
RPM_columns = list(filter(lambda col: col.startswith('RPM_'), df.columns))

df = df.astype({col_name: int for col_name in treatment_columns})
df = df.astype({col_name: float for col_name in RPM_columns})

RNaseR_polyAminus_col = list(filter(lambda col: col in ['RNaseR', 'polyA_minus'], treatment_columns))[0]
RNaseR_polyAminus_RPM_col = list(filter(lambda col: col.startswith(f'RPM_{RNaseR_polyAminus_col}'), RPM_columns))[0]

totalRNA_RPM_col = list(filter(lambda col: col.startswith('RPM_totalRNA'), RPM_columns))[0]

df = df.assign(
    **{
        f'Not detected({RNaseR_polyAminus_col})': (df[RNaseR_polyAminus_col] == 0).apply(int),
        f'Not depleted({RNaseR_polyAminus_col})': (df[totalRNA_RPM_col] <= df[RNaseR_polyAminus_RPM_col]).apply(int)
    }
)

df.to_csv(snakemake.output[0], index=False, sep='\t')
