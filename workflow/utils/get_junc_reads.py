#! /usr/bin/env python

import sys
import argparse
from sam import SamFormat


def get_junc_reads(ref_len, cross_junc_threshold, map_len_threshold, similarity_threshold):
    for line in sys.stdin:
        sam_data = SamFormat(line)

        if not sam_data.is_header:
            if not sam_data.is_unmapped:
                pos = sam_data.pos
                Z3 = int(sam_data.optional_fields['Z3'].value)
                XS = float(sam_data.optional_fields['XS'].value)

                map_len = Z3 - sam_data.pos + 1

                if (map_len >= map_len_threshold) and \
                        (XS >= similarity_threshold) and \
                        (pos <= ref_len - cross_junc_threshold + 1) and \
                        (Z3 >= ref_len + cross_junc_threshold):

                    print(
                        sam_data.qname,
                        sam_data.rname,
                        sam_data.pos,
                        Z3,
                        map_len,
                        sam_data.cigar,
                        sam_data.optional_fields['MD'].value,
                        XS,
                        sep='\t'
                    )


def create_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('ref_len', type=float)
    parser.add_argument('junction_dist', type=float)
    parser.add_argument('map_len', type=int)
    parser.add_argument('similarity', type=float)

    return parser


if __name__ == "__main__":
    parser = create_parser()
    args = parser.parse_args()

    get_junc_reads(args.ref_len, args.junction_dist, args.map_len, args.similarity)
