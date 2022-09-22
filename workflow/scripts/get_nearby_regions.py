import csv


with open(snakemake.input[0], newline='', encoding='utf-8') as f_in, \
        open(snakemake.output[0], 'w') as out:

    _ = f_in.readline()
    reader = csv.reader(f_in, delimiter='\t')

    for circRNA in reader:
        chr_ = f"chr{circRNA[0]}"
        pos1, pos2 = sorted(map(int, circRNA[1:3]))
        strand = circRNA[3]
        circ_id = circRNA[4]

        nearby_regions = [
            (pos1 - 10, pos1 - 1),
            (pos1, pos1 + 10 - 1),
            (pos2 - 10 + 1, pos2),
            (pos2 + 1, pos2 + 10)
        ]

        if strand == "+":
            print(chr_, *nearby_regions[2], circ_id, 'donor_in', sep='\t', file=out)
            print(chr_, *nearby_regions[3], circ_id, 'donor_out', sep='\t', file=out)
            print(chr_, *nearby_regions[1], circ_id, 'acceptor_in', sep='\t', file=out)
            print(chr_, *nearby_regions[0], circ_id, 'acceptor_out', sep='\t', file=out)
        else:
            print(chr_, *nearby_regions[1], circ_id, 'donor_in', sep='\t', file=out)
            print(chr_, *nearby_regions[0], circ_id, 'donor_out', sep='\t', file=out)
            print(chr_, *nearby_regions[2], circ_id, 'acceptor_in', sep='\t', file=out)
            print(chr_, *nearby_regions[3], circ_id, 'acceptor_out', sep='\t', file=out)
