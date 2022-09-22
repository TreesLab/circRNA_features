import argparse
import os
import os.path
import multiprocessing as mp
import time
import re

from selenium import webdriver
from RBPmap import RBPmapJob

import logging

logging.basicConfig(
    encoding='utf-8',
    format="{asctime} ({name}) [{levelname}] {message}",
    style='{',
    level=logging.INFO
)

def submit_jobs(jobs, out_files, queue):
    for job, out_file in zip(jobs, out_files):
        job.init_page()
        job.input_data()
        job.submit()

        logging.info(f'Job {job.job_id} was submitted! ({job.job_name})')

        job.init_page()

        queue.put((job.job_id, out_file))


def download_results(queue, total):
    results_count = 0

    while results_count < total:
        job_id, out_file = queue.get()

        logging.info(f'Waiting for job {job_id} to be finished.')

        while RBPmapJob.get_job_status(job_id) != 'Finished':
            time.sleep(60)

        logging.info(f'Job {job_id} has been finished.')

        RBPmapJob._save_result(job_id, out_file)
        logging.info(f'The result of job {job_id} was downloaded!')

        results_count += 1


def create_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_dir')
    parser.add_argument('output_dir')
    parser.add_argument('--prefix')
    parser.add_argument('--job_limit', type=int, default=5)
    parser.add_argument('--do-not-overwrite', action='store_true')
    parser.add_argument('--high-cons', action='store_true', help='Use "high stringency" and "conservation filter".')

    return parser


if __name__ == "__main__":
    parser = create_parser()
    args = parser.parse_args()

    mp.set_start_method('fork')

    os.makedirs(args.output_dir, exist_ok=True)

    in_files = list(
        filter(
            lambda entry: entry.name.startswith(args.prefix),
            os.scandir(args.input_dir)
        )
    )
    in_files = sorted(
        in_files,
        key=lambda entry: int(re.search(rf'(?<={args.prefix})[0-9]+', entry.name).group())
    )
    logging.debug(f'total {len(in_files)} input files.')

    logging.debug('generating the paths of output files...')
    out_file_paths = list(map(
        lambda entry: os.path.join(args.output_dir, f'{entry.name}.RBPmap_result'),
        in_files
    ))
    logging.debug('generating the paths of output files... done')


    # start driver
    driver = webdriver.Firefox()
    logging.info('The driver started!')

    all_RBPmap_jobs = [
        RBPmapJob(
            driver,
            entry.name,
            entry.path,
            high_cons=args.high_cons
        )
        for entry in in_files
    ]

    if args.do_not_overwrite:
        logging.info('user specifies "--do-not-overwrite"!')
        all_RBPmap_jobs_tmp = []
        out_file_paths_tmp = []

        for job, out_file in zip(all_RBPmap_jobs, out_file_paths):
            if not os.path.exists(out_file):
                all_RBPmap_jobs_tmp.append(job)
                out_file_paths_tmp.append(out_file)
            else:
                logging.info(f'The file {out_file} exists! Omitting the corresponding job.')

        all_RBPmap_jobs = all_RBPmap_jobs_tmp
        out_file_paths = out_file_paths_tmp

    logging.info(f'There are {len(all_RBPmap_jobs)} job(s) to be submitted.')

    job_queue = mp.Queue(args.job_limit)

    submitter = mp.Process(target=submit_jobs, args=(all_RBPmap_jobs, out_file_paths, job_queue))
    downloader = mp.Process(target=download_results, args=(job_queue, len(all_RBPmap_jobs)))


    submitter.start()
    time.sleep(30)
    downloader.start()

    submitter.join()
    logging.info('submitter finishes all jobs.')

    downloader.join()
    logging.info('downloader finishes all jobs.')

    logging.info('closing driver...')
    driver.close()
    logging.info('done.')
