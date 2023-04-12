"""
Input format:

1:633567|634095(-)      1       DDX3X   1.0     1
1:633567|634095(-)      1       AGO1-4  0.8466666666666667      239
1:633567|634095(-)      1       AGO1-4  0.8470588235294118      569
1:633567|634095(-)      1       AGO1-4  0.896   329
1:633567|634095(-)      1       DDX3X   1.0     1
1:633567|634095(-)      2       AGO1-4  0.823   1
1:633567|634095(-)      2       SRSF3   0.9333333333333333      1
1:633567|634095(-)      2       NOP56   1.0     4
1:633567|634095(-)      2       BUD13   1.0     13
1:633567|634095(-)      2       AGO1-4  1.0     14

"""

import csv
from itertools import groupby
from operator import itemgetter
from collections import defaultdict


def coverage_filter(data_reader, threshold=0.8):
    for data in data_reader:
        if float(data[3]) >= threshold:
            yield data


with open(snakemake.input[0]) as f_in, open(snakemake.output[0], 'w') as out:
    reader = coverage_filter(csv.reader(f_in, delimiter='\t'))

    for circ_id, circ_gp in groupby(reader, key=itemgetter(0)):
        RBP_data = {}
        for region_id, region_gp in groupby(circ_gp, key=itemgetter(1)):
            region_RBPs = defaultdict(lambda: 100000)
            for RBP, dist in map(itemgetter(2, 4), region_gp):
                if int(dist) < region_RBPs[RBP]:
                    region_RBPs[RBP] = int(dist)

            RBP_data[region_id] = region_RBPs

        if len(RBP_data) == 2:
            RBPs_1 = set(RBP_data['1'].keys())
            RBPs_2 = set(RBP_data['2'].keys())
            common_RBP = sorted(RBPs_1 & RBPs_2)

            if common_RBP:
                for RBP in common_RBP:
                    min_dist = min([RBP_data['1'][RBP], RBP_data['2'][RBP]])
                    print(circ_id, RBP, min_dist, sep='\t', file=out)
