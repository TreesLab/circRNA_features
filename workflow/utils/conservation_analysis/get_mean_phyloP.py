#! /usr/bin/env python

import argparse
import pyBigWig
import statistics as stat


def create_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-bw', dest='bigWig_file', help='phyloP bigWig file.',
                        required=True)
    parser.add_argument('regions_file', type=argparse.FileType('r'))

    return parser


if __name__ == "__main__":
    parser = create_parser()
    args = parser.parse_args()

    bw = pyBigWig.open(args.bigWig_file)

    for line in args.regions_file:
        data = line.rstrip('\n').split('\t')
        chr_ = data[0]
        pos1 = int(data[1])
        pos2 = int(data[2])

        phyloP_scores = bw.values(chr_, pos1 - 1, pos2)
        mean_phyloP_score = stat.fmean(phyloP_scores)

        print(*data, mean_phyloP_score, sep='\t')
