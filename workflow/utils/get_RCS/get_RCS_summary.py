#! /usr/bin/env python

import argparse
from itertools import groupby
from collections import namedtuple, Counter


INPUT_TITLES = [
    'circRNA_id',
    'region1',
    'region2',
    'matches',
    'alignment_len',
    'mismatch',
    'gapopen',
    'alignment_start_1',
    'alignment_end_1',
    'alignment_start_2',
    'alignment_end_2',
    'E_value',
    'bit_score',
    'region1_type',
    'region2_type'
]


OUTPUT_TITLES = [
    'circRNA_id',
    '#RCS_across_flanking_introns',
    'min_dist_of_RCS(across, donor)',
    'min_dist_of_RCS(across, acceptor)',
    'sum_of_the_min_dist',
    '#RCS_within_donor\'s_intron',
    '#RCS_within_acceptor\'s_intron'
]


class RCSReader:
    RCS_data = namedtuple('RCS_data', INPUT_TITLES)

    def __init__(self, RCS_result_file):
        self._RCS_result_file = RCS_result_file
        self._reader = self._RCS_data_reader()

    def _RCS_data_reader(self):
        self._title = self._RCS_result_file.readline()

        for line in self._RCS_result_file:
            data = line.rstrip('\n').split('\t')
            data = self.RCS_data(*data)
            yield data

    def __iter__(self):
        return self._reader


def create_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('RCS_file', type=argparse.FileType('r'))
    parser.add_argument('--dist', default=2000, type=int)

    return parser


if __name__ == "__main__":
    parser = create_parser()
    args = parser.parse_args()

    print(*OUTPUT_TITLES, sep='\t')

    reader = RCSReader(args.RCS_file)
    group_reader = groupby(reader, key=lambda data: data.circRNA_id)

    for circRNA_id, RCS_gp in group_reader:
        RCS_gp = list(RCS_gp)
        RCS_types = Counter(
            [(data.region1_type, data.region2_type) for data in RCS_gp]
        )

        # dist
        across_RCS_gp = list(
            filter(
                lambda data: (data.region1_type == 'donor') and (data.region2_type == 'acceptor'),
                RCS_gp
            )
        )

        if len(across_RCS_gp) > 0:
            min_dist_donor = min([int(data.alignment_start_1) for data in across_RCS_gp])
            min_dist_acceptor = min([args.dist - int(data.alignment_start_2) for data in across_RCS_gp])
        else:
            min_dist_donor = args.dist
            min_dist_acceptor = args.dist

        total_min_dist = min_dist_donor + min_dist_acceptor

        print(
            circRNA_id,
            RCS_types[('donor', 'acceptor')],
            min_dist_donor,
            min_dist_acceptor,
            total_min_dist,
            RCS_types[('donor', 'donor')],
            RCS_types[('acceptor', 'acceptor')],
            sep='\t'
        )
