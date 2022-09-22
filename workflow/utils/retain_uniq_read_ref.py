#! /usr/bin/env python

import sys
from itertools import groupby
from collections import namedtuple


junc_reads_data = namedtuple(
    'JuncReads',
    [
        'read_id',
        'ref_id',
        'map_start',
        'map_end',
        'map_len',
        'CIGAR',
        'MD',
        'similarity',
        'read_end'
    ]
)


def get_uniq_read_ref(data_gp):
    return list(set([(data.read_id, data.ref_id) for data in data_gp]))


def retain_uniq_read_ref():
    data_reader = map(
        lambda line: junc_reads_data(*line.rstrip('\n').split('\t')),
        sys.stdin
    )

    for k, gp in groupby(data_reader, key=lambda data: data.read_id):
        gp = list(gp)
        uniq_read_ref_data = get_uniq_read_ref(gp)

        if len(uniq_read_ref_data) == 1:
            print(*uniq_read_ref_data[0], sep='\t')
        else:
            pass


if __name__ == "__main__":
    retain_uniq_read_ref()
