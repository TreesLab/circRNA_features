import re
import subprocess as sp
from collections import namedtuple

import logging

logging.basicConfig(
    filename=snakemake.log[0],
    level=logging.INFO,
    style='{',
    format='{asctime} - {message}'
)


def parse_fasta(fasta_txt):
    for m in re.finditer(r">(.+?)\n([^>]+)", fasta_txt):
        fa_id = m.group(1)
        fa_seq = ''.join(m.group(2).rstrip('\n').split('\n'))
        yield [fa_id, fa_seq]


def qgrs(seq):
    cmd = ['qgrs', '-g', '@', '-notitle']
    with sp.Popen(cmd, stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.PIPE, encoding='utf-8') as p:
        stdout, stderr = p.communicate(seq)

    return stdout, stderr


qgrs_result = namedtuple('qgrsResult', ['ID', 'T1', 'T2', 'T3', 'T4', 'TS', 'GS', 'SEQ'])


with open(snakemake.input[0]) as f_in, open(snakemake.output[0], 'w') as out:
    print('circ_id', *qgrs_result._fields, sep='\t', file=out)

    fasta_txt = f_in.read()

    for fa_id, fa_seq in parse_fasta(fasta_txt):
        results, err_msg = qgrs(fa_seq)
        logging.debug(fa_id)
        logging.debug(results)

        if not err_msg:
            for result in results.strip('\n').split('\n'):
                result = re.split(r' +', result)
                logging.debug(result)

                result = qgrs_result._make(result)
                logging.debug(result)

                print(fa_id, *result, sep='\t', file=out)
