import csv


with open(snakemake.input[0], newline='', encoding='utf-8') as f_in, \
        open(snakemake.output[0], 'w') as out:

    reader = csv.reader(f_in, delimiter='\t')

    for data in reader:
        L1 = int(data[2]) - int(data[1])
        L2 = int(data[8]) - int(data[7])
        L = min(L1, L2)

        if L > 0:
            coverage = float(data[12]) / L
        else:
            coverage = 0

        print(*data, L1, L2, L, coverage, sep='\t', file=out)
