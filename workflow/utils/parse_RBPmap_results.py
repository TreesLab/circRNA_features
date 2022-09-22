#! /usr/bin/env python

import argparse
import re


def parse_raw_result(raw_RBPmap_result_text):
    results = re.finditer(
        r'(?P<seq_id>[^\n]+)\n=+\n(?P<all_rbp_results>[^=]+?)\*+',
        raw_RBPmap_result_text,
        flags=re.DOTALL
    )

    for m in results:
        yield m.group('seq_id'), m.group('all_rbp_results')


_aln_pat = re.compile(
    r'^([^ \t]+)(?:[ \t]+([^ \t]+))?[ \t]+([^ \t]+)[ \t]+([^ \t]+)[ \t]+([^ \t]+)[ \t]+([^ \t]+)$',
    flags=re.M
)

def parse_rbp_results(rbp_results_text):
    rbp_results = re.split(r'\n(?=Protein:)', rbp_results_text.strip('\n'))
    for rbp_result in rbp_results:
        m = re.search(
            r'Protein: (?P<rbp>[^\n]+)\((?P<species>[^\n]+)\)\n(?P<all_aln>.+)',
            rbp_result.strip('\n'),
            flags=re.DOTALL
        )
        if m:
            all_aln_data = re.findall(_aln_pat, m.group('all_aln'))

            for data in all_aln_data[1:]:
                output = (m.group('rbp'), m.group('species'), *data)

                yield output


def create_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('RBPmap_result', nargs='+')

    return parser


if __name__ == "__main__":
    parser = create_parser()
    args = parser.parse_args()

    try:
        titles = [
            'seq_id',
            'rbp',
            'species',
            'position',
            'genomic_coordinate',
            'motif',
            'k_mer',
            'z_score',
            'p_value'
        ]

        print(*titles, sep='\t')

        for file_ in args.RBPmap_result:
            with open(file_) as f_in:
                raw_RBPmap_result_text = f_in.read()

            results = parse_raw_result(raw_RBPmap_result_text)
            for seq_id, all_rbp_results_text in results:
                all_rbp_results = parse_rbp_results(all_rbp_results_text)
                for rbp_result in all_rbp_results:
                    print(seq_id, *rbp_result, sep='\t')

    except BrokenPipeError:
        pass
