import pandas as pd
import io
import re
from itertools import groupby
from operator import itemgetter

import logging

logging.basicConfig(
    filename=snakemake.log[0],
    encoding='utf-8',
    format="{asctime} ({name}) [{levelname}] {message}",
    style='{',
    level=logging.INFO
)


def parse_line_data(line):
    return line.rstrip('\n').split('\t')

L = snakemake.params.L
cross_junction = snakemake.params.cross_junction

with open(snakemake.input[0]) as f_in, \
        open(snakemake.output[0], 'w') as out_1, open(snakemake.output[1], 'w') as out_2:

    titles = parse_line_data(f_in.readline())
    logging.debug(titles)

    print(*titles, 'len_motif', 'pos_end', 'cross_junction', sep='\t', file=out_1)
    print('seq_id', 'count', 'min_p_value', 'RBPs', sep='\t', file=out_2)

    for seq_id, gp in groupby(f_in, key=lambda line: re.search(r'^([^\t]+)', line).group(1)):
        with io.StringIO() as temp:
            print(*titles, sep='\t', file=temp)
            print(*gp, sep='', end='', file=temp)
            temp.seek(0)

            logging.debug(temp.getvalue())

            df = pd.read_csv(temp, sep='\t', dtype='object')
            logging.debug(df.columns)

            df = df.astype({'position': int, 'p_value': float})

            df = df[df['p_value'] <= 0.005]

            df = df.assign(len_motif=df['motif'].apply(len))
            df = df.assign(pos_end=df['position'] + df['len_motif'])
            df = df.assign(
                cross_junction=((df['position'] <= L - cross_junction + 1) & (df['pos_end'] >= L + cross_junction)).apply(int)
            )

            print(df.to_csv(sep='\t', index=False, header=False), end='', file=out_1)

            count_df = df[df['cross_junction'] == 1].groupby(
                'seq_id'
            ).agg(
                {
                    'rbp': [
                        lambda rbps: len(set(rbps)),
                        lambda rbps: ','.join(sorted(set(rbps)))
                    ],
                    'p_value': 'min'
                }
            ).reset_index()
            count_df.columns = ['seq_id', 'count', 'RBPs', 'min_p_value']
            count_df = count_df[['seq_id', 'count', 'min_p_value', 'RBPs']]

            print(count_df.to_csv(sep='\t', index=False, header=False), end='', file=out_2)
