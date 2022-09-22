#! /usr/bin/env python

"""
Usage:
    cat data.tsv | python ./merge_rows.py > out.tsv


Input:

K1    A
K2    B
K2    C
K3    D

--

Output:

K1    1    A
K2    2    B,C
K3    1    D

"""

import sys
from itertools import groupby
from operator import itemgetter


def data_reader(file_obj):
    for line in file_obj:
        data = line.rstrip('\n').split('\t')
        yield data


if __name__ == "__main__":
    data = data_reader(sys.stdin)
    for k, gp in groupby(data, key=itemgetter(0)):
        gp_data = list(map(itemgetter(1), gp))
        print(k, len(gp_data), ','.join(gp_data), sep='\t')
