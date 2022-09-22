import requests
import re
import tempfile as tp
import os.path
import time

import logging

logging.basicConfig(
    filename=snakemake.log[0],
    level=logging.INFO,
    format='{asctime}:{levelname}:{message}',
    style='{'
)
logging.info('test')

def read_fasta_file(fasta_file):
    with open(fasta_file) as f_in:
        fasta_data = re.findall(r'^>[^>]+$', f_in.read(), flags=re.M)
        fasta_data = list(map(lambda line: line.rstrip('\n'), fasta_data))

    return fasta_data


def save_fasta_to_file(fasta_data, fasta_file):
    logging.info('saving fasta data to file...')
    with open(fasta_file, 'w') as out:
        print('\n'.join(fasta_data), file=out)
    logging.info('done!')


def get_chunk(data_list, chunk_size):
    for i in range(0, len(data_list), chunk_size):
        yield data_list[i:i + chunk_size]


def request_MaxEntScan_5p(fasta_file):
    logging.info('make request to MaxEntScan_5p...')
    r = requests.post(
        'http://hollywood.mit.edu/cgi-bin/Xmaxentscan_scoreseq.pl',
        data={
            'MAXENT': 1,
            'MM': 1,
            'WMM': 1,
        },
        files={
            'INPUTFILE': (
                os.path.basename(fasta_file),
                open(fasta_file, 'rb')
            )
        }
    )
    logging.info(f'response status: {r.status_code}')

    return r


def request_MaxEntScan_3p(fasta_file):
    logging.info('make request to MaxEntScan_3p...')
    r = requests.post(
        'http://hollywood.mit.edu/cgi-bin/Xmaxentscan_scoreseq_acc.pl',
        data={
            'MAXENT': 1,
            'MM': 1,
            'WMM': 1,
        },
        files={
            'INPUTFILE': (
                os.path.basename(fasta_file),
                open(fasta_file, 'rb')
            )
        }
    )
    logging.info(f'response status: {r.status_code}')

    return r


def parse_data_from_MaxEntScan(MaxEntScan_output):
    data = re.findall(
        r'^>(.+)\([+-]\)\n([ATCG]+).*MAXENT: ([-.0-9]+).*MM: ([-.0-9]+).*WMM: ([-.0-9]+).*',
        MaxEntScan_output,
        flags=re.M
    )

    if data:
        return data
    else:
        raise Exception("Empty results!")


# five prime
logging.info('Five prime part.')

five_prime_fasta = read_fasta_file(snakemake.input.five_prime_fa)

with open(snakemake.output.five_prime_results, 'w') as out:
    with tp.NamedTemporaryFile(dir='.', prefix="MaxEntScan_tmp.", suffix='.fa') as tmp_file, \
            tp.TemporaryDirectory() as tmp_dir:
        for i, chunk in enumerate(get_chunk(five_prime_fasta, 10000), start=1):
            logging.info(f'---------- Chunk {i} ----------')
            save_fasta_to_file(chunk, tmp_file.name)
            r = request_MaxEntScan_5p(tmp_file.name)
            results = parse_data_from_MaxEntScan(r.text)
            logging.info(f'#results = {len(results)}')

            logging.info(f'output results of chunk {i}')
            for res in results:
                print(*res, sep='\t', file=out)

            time.sleep(0.05)


# three prime
logging.info('Three prime part')

three_prime_fasta = read_fasta_file(snakemake.input.three_prime_fa)

with open(snakemake.output.three_prime_results, 'w') as out:
    with tp.NamedTemporaryFile(dir='.', prefix="MaxEntScan_tmp.", suffix='.fa') as tmp_file:
        for i, chunk in enumerate(get_chunk(three_prime_fasta, 10000), start=1):
            logging.info(f'---------- Chunk {i} ----------')
            save_fasta_to_file(chunk, tmp_file.name)
            r = request_MaxEntScan_3p(tmp_file.name)
            results = parse_data_from_MaxEntScan(r.text)
            logging.info(f'#results = {len(results)}')

            logging.info(f'output results of chunk {i}')
            for res in results:
                print(*res, sep='\t', file=out)

            time.sleep(0.05)
