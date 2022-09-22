#! /usr/bin/env python


import argparse
import multiprocessing as mp
import subprocess as sp
import os.path


class CheckAAReads:
    def __init__(self, ref_genome, ref_others, checkAA_bin='checkAA_reads.py'):
        self.ref_genome = ref_genome
        self.ref_others = ref_others
        self.checkAA_bin = checkAA_bin

    def run(self, in_file, out_file):
        cmd = self._generate_cmd(in_file, out_file)
        res = sp.run(cmd)
        return res

    def _generate_cmd(self, in_file, out_file):
        cmd = [self.checkAA_bin]
        cmd += ['-rg', self.ref_genome]
        cmd += ['-ro', self.ref_others]
        cmd += [in_file, out_file]

        return cmd


def create_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-rg', dest='ref_genome', required=True)
    parser.add_argument('-ro', dest='ref_others', required=True)
    parser.add_argument('file_list', type=argparse.FileType('r'))
    parser.add_argument('out_dir')
    parser.add_argument('-p', '--num_proc', type=int, default=1)
    parser.add_argument('--checkAA_bin', default='checkAA_reads.py')

    return parser


if __name__ == "__main__":
    parser = create_parser()
    args = parser.parse_args()

    file_list = [filename.rstrip('\n') for filename in args.file_list]
    out_file_list = [
        os.path.join(
            args.out_dir,
            f'{os.path.basename(filename)}.checkAA_result'
        )
        for filename in file_list
    ]

    checkAA = CheckAAReads(
        args.ref_genome,
        args.ref_others,
        args.checkAA_bin
    )

    with mp.Pool(processes=args.num_proc) as pool:
        checkAA_results = [
            pool.apply_async(
                checkAA.run,
                (in_file, out_file)
            )
            for in_file, out_file in zip(file_list, out_file_list)
        ]

        for res in checkAA_results:
            res.wait()
