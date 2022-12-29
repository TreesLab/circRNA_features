import csv
import logging

logging.basicConfig(
    filename=snakemake.log[0],
    level=logging.DEBUG,
    style='{',
    format='{asctime} - {message}'
)


def is_across_junction(data, junc_pos, dist):
    structure_start = int(data['T1'])
    structure_end = int(data['T4']) + int(data['TS'])

    if structure_start <= junc_pos - dist + 1:
        if structure_end >= junc_pos + dist:
            return True

    return False


with open(snakemake.input[0]) as f_in, open(snakemake.output[0], 'w') as out:
    title = next(f_in)
    title = title.rstrip('\n').split('\t')
    print(*title, sep='\t', file=out)
    logging.debug(title)

    reader = csv.DictReader(f_in, delimiter='\t', fieldnames=title)

    for data in reader:
        logging.debug(data)

        if is_across_junction(data, snakemake.params.L - 1, int(snakemake.wildcards.dist)):
            print(*data.values(), sep='\t', file=out)
