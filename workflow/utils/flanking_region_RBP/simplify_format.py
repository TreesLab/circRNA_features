#!    /usr/bin/env    python

"""
# Usage:

cat human.circRNAs.totalRNA.tsv.flanking_1k.bed.intersect_RBP.raw.coverage | python simplify_format.py > out.tsv


# Example

File:
human.circRNAs.totalRNA.tsv.flanking_1k.bed.intersect_RBP.raw.coverage

Input:
1    632566    633566    1:633567|634095(-)_1    .    -    1    631678    632668    SBDH81_AGO1-4_13    255    +    102    1000    990    990    0.10303
1    632566    633566    1:633567|634095(-)_1    .    -    1    631737    634372    SBDH439_DDX3X_2    255    -    1000    1000    2635    1000    1
1    632566    633566    1:633567|634095(-)_1    .    -    1    631738    633328    SBDH82_AGO1-4_7    255    +    762    1000    1590    1000    0.762
1    634095    635095    1:633567|634095(-)_2    .    -    1    631737    634372    SBDH439_DDX3X_2    255    -    277    1000    2635    1000    0.277
1    634095    635095    1:633567|634095(-)_2    .    -    1    632521    634322    SBDH438_DDX3X_2    255    -    227    1000    1801    1000    0.227
1    634095    635095    1:633567|634095(-)_2    .    -    1    633358    634438    SBDH76_AGO1-4_19    255    +    343    1000    1080    1000    0.343

Output:
1:633567|634095(-)    1    AGO1-4    0.10303
1:633567|634095(-)    1    DDX3X     1
1:633567|634095(-)    1    AGO1-4    0.762
1:633567|634095(-)    2    DDX3X    0.277
1:633567|634095(-)    2    DDX3X    0.227
1:633567|634095(-)    2    AGO1-4    0.343
"""

import sys
import re


if __name__ == "__main__":
    for line in sys.stdin:
        data = line.rstrip('\n').split('\t')

        region_id, pair_id = re.split(r'_(?=[12]$)', data[3])

        if data[9] == '.':
            RBP_name = '.'
        else:
            RBP_name = data[9].split('_')[1]

        coverage = data[16]

        print(region_id, pair_id, RBP_name, coverage, sep='\t')
