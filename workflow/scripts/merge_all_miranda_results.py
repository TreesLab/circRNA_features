import os
import io
import re
import logging

logging.basicConfig(filename=snakemake.log[0], encoding='utf-8', level=logging.DEBUG)


in_files = list(
    filter(
        lambda entry: entry.name.startswith(snakemake.params.filename_prefix),
        os.scandir(snakemake.input[0])
    )
)
logging.info(f'There are {len(in_files)} in_files.')
logging.debug("\n".join([file_.name for file_ in in_files[:10]]))

in_files = sorted(
    in_files,
    key=lambda entry: int(re.search(rf'(?<={snakemake.params.filename_prefix})[0-9]+', entry.name).group())
)
logging.debug("\n".join([file_.name for file_ in in_files[:10]]))


RESULT_TITLE = (
    'query_id',
    'reference_id',
    'score',
    'energy',
    'query_start',
    'query_end',
    'ref_start',
    'ref_end',
    'aln_length',
    'identity',
    'similarity',
    'aln_mirna',
    'aln_map',
    'aln_utr'
)

logging.info('Merging results ...')
with open(snakemake.output[0], 'w') as out:
    print(*RESULT_TITLE, sep='\t', file=out)

    for file_ in in_files:
        with open(file_.path) as f_in:
            print(f_in.read(), file=out)

logging.info('Done.')