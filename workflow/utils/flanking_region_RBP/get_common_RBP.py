#! /usr/bin/env python

"""
Input:
1:633567|634095(-)    1    AGO1-4    0.10303
1:633567|634095(-)    1    DDX3X     1
1:633567|634095(-)    1    AGO1-4    0.762
1:633567|634095(-)    2    DDX3X    0.277
1:633567|634095(-)    2    DDX3X    0.227
1:633567|634095(-)    2    AGO1-4    0.343

Output:
1:633567|634095(-)    2    AGO1-4,DDX3X
"""

import sys
from itertools import groupby
from collections import namedtuple


Data = namedtuple('Data', ('region_id', 'pair_id', 'RBP', 'coverage'))


def get_data():
    for line in sys.stdin:
        data = line.rstrip('\n').split('\t')
        data[3] = float(data[3])
        data = Data(*data)

        yield data


if __name__ == "__main__":
    for region_id, gp in groupby(get_data(), key=lambda data: data.region_id):
        gp = list(filter(lambda data: data.coverage >= 0.8, gp))
        RBP_1 = set(map(lambda data: data.RBP, filter(lambda data: data.pair_id == '1', gp)))
        RBP_2 = set(map(lambda data: data.RBP, filter(lambda data: data.pair_id == '2', gp)))
        common_RBP = tuple(sorted(RBP_1 & RBP_2))

        print(region_id, len(common_RBP), ','.join(common_RBP), sep='\t')
