import os
import os.path
import re
import multiprocessing as mp
import logging

logging.basicConfig(filename=snakemake.log[0], encoding='utf-8', level=logging.INFO)


in_files = list(
    filter(
        lambda entry: entry.name.startswith(snakemake.params.filename_prefix),
        os.scandir(snakemake.input[0])
    )
)
logging.info(f'There are {len(in_files)} in_files.')
logging.debug(in_files)


os.makedirs(snakemake.output[0])
logging.info('Output directory created successfully!')


out_file_paths = list(map(
    lambda entry: os.path.join(snakemake.output[0], f'{entry.name}.parsed'),
    in_files
))
logging.debug(out_file_paths)


def _get_values(res_line):
    values = [kv.split('=')[1] for kv in res_line.split('\t')]
    logging.debug(values)
    return values

def _parse_result(raw_result):
    aln_res_pat = re.compile(r'^//hit_info\t(.*)', flags=re.M)
    for m in re.finditer(aln_res_pat, raw_result):
        logging.debug(m[1])
        yield _get_values(m.group(1))


def parse_miranda_results(in_file_path, out_file_path):
    with open(in_file_path) as f_in, \
            open(out_file_path, 'w') as out:

        raw_result_text = f_in.read()

        logging.debug(raw_result_text[:5000])

        for res in _parse_result(raw_result_text):
            print(*res, sep='\t', file=out)
            logging.debug(res)

logging.info('Start to parse the miranda results.')
with mp.Pool(processes=snakemake.threads) as pool:
    processes = [
        pool.apply_async(
            parse_miranda_results,
            (in_file.path, out_file_path)
        )
        for in_file, out_file_path in zip(in_files, out_file_paths)
    ]
    logging.info('All processes created!')

    logging.info('Waiting ...')

    for p in processes:
        logging.debug(f'Process {p} completed!')
        p.wait()

    logging.info('All processes completed!')
