#! /usr/bin/env python

"""
Example.
    python merge_junc_reads_result.py \
        human.circRNAs.tsv \
        totalRNA/all_junc_reads.count \
        RNaseR/all_junc_reads.count \
        -c TotalRNA RNaseR \
        > out.tsv

"""


import argparse
import pandas as pd


def merge_result(circRNAs_file, data_files, out_file=None, column_names=None):
    circ_df = pd.read_csv(circRNAs_file, sep='\t', dtype='object')

    if column_names is None:
        column_names = (f"data{i}" for i in range(1, len(data_files) + 1))

    for data_file, col_name in zip(data_files, column_names):
        data_df = pd.read_csv(
            data_file,
            sep='\t',
            names=['event_id', col_name],
            dtype='object',
            usecols=[0, 1]
        )

        circ_df = circ_df.merge(data_df, on='event_id', how='left').fillna({col_name: '0'})

    if out_file is None:
        print(circ_df.to_csv(sep='\t', index=False))
    else:
        circ_df.to_csv(out_file, sep='\t', index=False)


def create_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('circRNAs')
    parser.add_argument('data_files', metavar='data_file', nargs='+')
    parser.add_argument('-c', '--column_names', nargs='+',
                        help='The column names of the data_files, resp.')
    parser.add_argument('-o', '--out_file')

    return parser


if __name__ == "__main__":
    parser = create_parser()
    args = parser.parse_args()

    merge_result(
        args.circRNAs,
        args.data_files,
        args.out_file,
        args.column_names
    )
