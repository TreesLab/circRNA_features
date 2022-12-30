"""
1:633567|634095(-)      ADAR    46
1:633567|634095(-)      AGO1    11
1:633567|634095(-)      AGO1-4  1
1:633567|634095(-)      AGO2    10
1:633567|634095(-)      AGO3    70
1:633567|634095(-)      ELAVL1  31
1:633567|634095(-)      FTO     1
1:633567|634095(-)      FUS     37
1:633567|634095(-)      GTF2F1  200
1:633567|634095(-)      HNRNPA1 354
"""

import csv
from itertools import groupby
from operator import itemgetter


with open(snakemake.input[0]) as f_in, open(snakemake.output[0], 'w') as out:
    reader = csv.reader(f_in, delimiter='\t')

    for circ_id, circ_gp in groupby(reader, key=itemgetter(0)):
        circ_gp = list(circ_gp)
        min_dist = min(int(data[2]) for data in circ_gp)
        RBP_list = ','.join(map(itemgetter(1), circ_gp))
        print(circ_id, min_dist, RBP_list, sep='\t', file=out)
