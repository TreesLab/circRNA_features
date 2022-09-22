#! /usr/bin/env python


import argparse
import os
import os.path
import re
import tempfile as tp
import shutil
from contextlib import contextmanager


def parse_fasta(fasta_txt):
    for m in re.finditer(r">(.+?)(?:\([+-]\))?\n([^>]+)", fasta_txt):
        fa_id = m.group(1)
        fa_seq = ''.join(m.group(2).rstrip('\n').split('\n'))
        yield [fa_id, fa_seq]


def generate_ensembl_protein_coding(cdna_file, out_file):
    with open(cdna_file) as fa_in:
        fa_txt = fa_in.read()

    with open(out_file, 'w') as out:
        for fa_id, fa_seq in parse_fasta(fa_txt):
            m = re.search(r'{}'.format("protein_coding"), fa_id)
            if m:
                print(">{}\n{}".format(fa_id, fa_seq), file=out)

    return out_file


def generate_ensembl_lncRNA(ncrna_file, out_file):
    with open(ncrna_file) as fa_in:
        fa_txt = fa_in.read()

    with open(out_file, 'w') as out:
        for fa_id, fa_seq in parse_fasta(fa_txt):
            m = re.search(r'{}'.format("lncRNA"), fa_id)
            if m:
                print(">{}\n{}".format(fa_id, fa_seq), file=out)

    return out_file


def generate_RepChrM(genome_file, out_file):
    with open(genome_file) as fa_in:
        fa_txt = fa_in.read()

    with open(out_file, 'w') as out:
        for fa_id, fa_seq in parse_fasta(fa_txt):
            m = re.search(r'^(?:chr)?MT?', fa_id)
            if m:
                print(">{}\n{}\n{}".format(fa_id, fa_seq, fa_seq), file=out)

    return out_file


@contextmanager
def cwd(path):
    origin_pwd = os.getcwd()
    os.chdir(path)
    yield
    os.chdir(origin_pwd)


def generate_other_ref(genome_file,
                       pc_src_file,
                       lncRNA_src_file,
                       out_file,
                       src_type='ensembl'):

    genome_file = os.path.abspath(genome_file)
    pc_src_file = os.path.abspath(pc_src_file)
    lncRNA_src_file = os.path.abspath(lncRNA_src_file)
    out_file = os.path.abspath(out_file)

    with tp.TemporaryDirectory(prefix='other_ref_tmp.', dir='.') as tmp_dir:
        with cwd(tmp_dir):
            if src_type == 'ensembl':
                pc_file = generate_ensembl_protein_coding(
                    pc_src_file,
                    'protein_coding.fa'
                )
                lncRNA_file = generate_ensembl_lncRNA(
                    lncRNA_src_file,
                    'lncRNA.fa'
                )
            elif src_type == 'gencode':
                pc_file = pc_src_file
                lncRNA_file = lncRNA_src_file

            repChrM_file = generate_RepChrM(genome_file, 'RepChrM.fa')

            with open(out_file, 'wb') as out:
                for ref_file in [pc_file, lncRNA_file, repChrM_file]:
                    with open(ref_file, 'rb') as f_in:
                        shutil.copyfileobj(f_in, out)


def create_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('genome_file')
    parser.add_argument('pc_src_file')
    parser.add_argument('lncRNA_src_file')
    parser.add_argument('out_file')
    parser.add_argument('--src_type', default='ensembl')

    return parser


if __name__ == "__main__":
    parser = create_parser()
    args = parser.parse_args()

    generate_other_ref(
        args.genome_file,
        args.pc_src_file,
        args.lncRNA_src_file,
        args.out_file,
        args.src_type
    )
