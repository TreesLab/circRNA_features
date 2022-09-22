#! /usr/bin/env python

import sys
import re
from sam import SamFormat


def append_XS_tag():
    for line in sys.stdin:
        sam_data = SamFormat(line)

        if sam_data.is_header:
            print(sam_data)
        else:
            if sam_data.is_unmapped:
                print(sam_data)
            else:
                XS_tag = generate_XS_tag(sam_data)
                print(sam_data, XS_tag, sep='\t')


def generate_XS_tag(sam_data):
    similarity = calc_similarity(
        sam_data.cigar,
        sam_data.optional_fields['MD'].value
    )
    return "XS:f:{}".format(round(similarity, 4))


def calc_similarity(cigar, MD):
    total_matches = sum(cigar['M'])
    perfect_matches = sum(map(int, re.findall(r'([0-9]+)', MD)))

    similarity = perfect_matches / total_matches

    return similarity


if __name__ == "__main__":
    append_XS_tag()
