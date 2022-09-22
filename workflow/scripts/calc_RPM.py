import pandas as pd

df = pd.read_csv(snakemake.input[0], sep='\t', dtype='object')

def get_num_reads(num_reads_file):
    with open(num_reads_file) as f_in:
        num_reads = int(f_in.read().rstrip('\n'))
        return num_reads

all_num_reads = [get_num_reads(file_) for file_ in snakemake.input[1:]]


column_names = [
    f'RPM_{name}(num_reads={num_reads})'
    for name, num_reads in zip(snakemake.params.treatments, all_num_reads)
]


df = df.assign(
    **{
        column_names[i]: df[snakemake.params.treatments[i]].astype(int) * 1000000 / all_num_reads[i]
        for i in range(len(column_names))
    }
)


df.to_csv(snakemake.output[0], index=False, sep='\t')
