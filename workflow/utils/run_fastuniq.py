#! /usr/bin/env python

"""
A python wrapper for fastuniq.

Usage:
  run_fastuniq.py [FQ_1] [FQ_2] [UNIQ_FQ_1] [UNIQ_FQ_2]

  FQ_1, FQ_2:

    input gzipped fastq files

  UNIQ_FQ_1, UNIQ_FQ_2:

    output fastq files with gzip

"""

# TODO
#   - [ ] use 'argparse' to parse arguments
#   - [ ] provide an option to retain the tmp_dir
#   - [v] add logging


import sys
import tempfile as tp
import subprocess as sp
import multiprocessing as mp
import re
import os
import os.path
import logging
from contextlib import contextmanager

logging.basicConfig(
    format="{asctime} - {message}",
    level=logging.INFO,
    style='{'
)


def print_usage():
    print(__doc__, file=sys.stderr)


def unzip_file(in_file, out_file):
    unzip_proc = sp.Popen(['zcat', in_file], stdout=sp.PIPE, encoding='utf-8')

    logging.info(f"Start to unzip the file: {in_file}")

    with open(out_file, 'w') as out:
        for line in unzip_proc.stdout:
            out.write(line)

    unzip_proc.wait()

    logging.info(f"The unzip process is completed! ({in_file})")


def unzip_files(in_files, out_files):
    unzip_procs = [
        mp.Process(target=unzip_file, args=(in_file, out_file))
        for in_file, out_file in zip(in_files, out_files)
    ]

    for proc in unzip_procs:
        proc.start()

    for proc in unzip_procs:
        proc.join()


def zip_file(in_file, out_file):
    zip_proc = sp.Popen(['gzip', '--fast', '-c', in_file], stdout=sp.PIPE)

    logging.info(f"Start to zip the file: {in_file}")

    with open(out_file, 'wb') as out:
        # out.write(zip_proc.stdout.read())
        for line in zip_proc.stdout:
            out.write(line)

    zip_proc.wait()

    logging.info(f"The zip process is completed! ({in_file})")


def zip_files(in_files, out_files):
    zip_procs = [
        mp.Process(target=zip_file, args=(in_file, out_file))
        for in_file, out_file in zip(in_files, out_files)
    ]

    for proc in zip_procs:
        proc.start()

    for proc in zip_procs:
        proc.join()


@contextmanager
def cwd(path):
    origin_pwd = os.getcwd()
    os.chdir(path)
    yield
    os.chdir(origin_pwd)


def fastuniq(fq_1, fq_2, uniq_fq_1, uniq_fq_2):
    logging.info(f"Starting the run_fastuniq process on {fq_1} & {fq_2}")

    fq_1 = os.path.abspath(fq_1)
    fq_2 = os.path.abspath(fq_2)
    uniq_fq_1 = os.path.abspath(uniq_fq_1)
    uniq_fq_2 = os.path.abspath(uniq_fq_2)

    with tp.TemporaryDirectory(prefix='fastuniq_tmp.', dir='.') as tmp_dir:
        with cwd(tmp_dir):
            # unzip files
            unzipped_fq_1 = re.sub(r'\.gz$', '', os.path.basename(fq_1))
            unzipped_fq_2 = re.sub(r'\.gz$', '', os.path.basename(fq_2))
            unzip_files([fq_1, fq_2], [unzipped_fq_1, unzipped_fq_2])

            # fastuniq
            fq_list = 'fq_list.txt'
            with open(fq_list, 'w') as out:
                print(unzipped_fq_1, file=out)
                print(unzipped_fq_2, file=out)

            unzipped_uniq_fq_1 = re.sub(r'\.(?:fastq|fq)$', '.uniq.fq', unzipped_fq_1)
            unzipped_uniq_fq_2 = re.sub(r'\.(?:fastq|fq)$', '.uniq.fq', unzipped_fq_2)

            logging.info(f"Start to run the fastuniq...")
            sp.run(['fastuniq', '-i', fq_list, '-o', unzipped_uniq_fq_1, '-p', unzipped_uniq_fq_2], stderr=sys.stderr)
            logging.info(f"The fastuniq process is completed!")

            # zip files
            zip_files([unzipped_uniq_fq_1, unzipped_uniq_fq_2], [uniq_fq_1, uniq_fq_2])

    logging.info(f"================ Process Completed ================")


if __name__ == "__main__":

    if len(sys.argv) < 5:
        print_usage()
        exit(1)

    fastuniq(*sys.argv[1:])
