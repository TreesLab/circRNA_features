import subprocess as sp
import multiprocessing as mp
import os
import os.path
import logging

logging.basicConfig(filename=snakemake.log[0], encoding='utf-8', level=logging.DEBUG)


os.makedirs(snakemake.output[0])


def run_miranda(ref_file, query_file, out_file, miranda_bin='miranda'):
    cmd = [miranda_bin, ref_file, query_file]
    cmd += ['-out', out_file]
    cmd += ['-keyval', '-quiet']

    res = sp.run(cmd, stderr=sp.PIPE, encoding='utf-8')

    logging.debug(res.stderr)

    return res


filenames = list(
    filter(
        lambda name: name.startswith(snakemake.params.filename_prefix),
        os.listdir(snakemake.input[1])
    )
)

in_files = list(map(
    lambda name: os.path.join(snakemake.input[1], name),
    filenames
))

out_files = list(map(
    lambda name: os.path.join(snakemake.output[0], f'{name}.miranda_result'),
    filenames
))


with mp.Pool(processes=snakemake.threads) as pool:
    miranda_results = [
        pool.apply_async(
            run_miranda,
            (snakemake.input[0], in_file, out_file)
        )
        for in_file, out_file in zip(in_files, out_files)
    ]

    for res in miranda_results:
        res.wait()
