import csv

dist = snakemake.params.dist

with open(snakemake.input[0], newline='', encoding='utf-8') as f_in, \
        open(snakemake.output[0], 'w') as out:

    _ = f_in.readline()
    reader = csv.reader(f_in, delimiter='\t')

    for circRNA in reader:
        pos = sorted(map(int, circRNA[1:3]))

        print(
            circRNA[0],
            pos[0] - 1 - dist,
            pos[0] - 1,
            f"{circRNA[4]}_1",
            0,
            circRNA[3],
            sep='\t',
            file=out
        )
        print(
            circRNA[0],
            pos[1] + 1 - 1,
            pos[1] + 1 - 1 + dist,
            f"{circRNA[4]}_2",
            0,
            circRNA[3],
            sep='\t',
            file=out
        )
