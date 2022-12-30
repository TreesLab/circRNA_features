"""
Input format:

1   632566  633566  1:633567|634095(-)_1    0   -   1   631737  634372  SBDH439_DDX3X_2 255 -   1000    1000    2635    1000    1.0
1   632566  633566  1:633567|634095(-)_1    0   -   1   632428  633328  SBDH79_AGO1-4_12    255 +   762 1000    900 900 0.8466666666666667
1   632566  633566  1:633567|634095(-)_1    0   -   1   632488  632998  SBDH77_AGO1-4_28    255 +   432 1000    510 510 0.8470588235294118
1   632566  633566  1:633567|634095(-)_1    0   -   1   632488  633238  SBDH80_AGO1-4_13    255 +   672 1000    750 750 0.896
1   632566  633566  1:633567|634095(-)_1    0   -   1   632521  634322  SBDH438_DDX3X_2 255 -   1000    1000    1801    1000    1.0

"""

class FlankingRegion:
    def __init__(self, flanking_region_data):
        self._parse_data(flanking_region_data)

    def _parse_data(self, data):
        self.chrm = data[0]
        self.start = int(data[1]) + 1
        self.end = int(data[2])
        self.circ_id, self.region_id = data[3].split('_')
        self.strand = data[5]

    @property
    def side(self):
        if self.region_id == '1':
            return 'LHS'
        elif self.region_id == '2':
            return 'RHS'


class RBPdata:
    def __init__(self, rbp_data):
        self._parse_data(rbp_data)

    def _parse_data(self, data):
        self.chrm = data[0]
        self.start = int(data[1]) + 1
        self.end = int(data[2])
        self.exp_id, self.RBP_name, _ = data[3].split('_')
        self.strand = data[5]



with open(snakemake.input[0]) as f_in, open(snakemake.output[0], 'w') as out:
    for line in f_in:
        data = line.rstrip('\n').split('\t')

        region = FlankingRegion(data[:6])
        coverage = float(data[16])

        if data[6] != '.':
            RBP_data = RBPdata(data[6:12])

            if region.side == 'LHS':
                RBP_dist = max(region.end - RBP_data.end + 1, 1)
            elif region.side == 'RHS':
                RBP_dist = max(RBP_data.start - region.start + 1, 1)

            print(region.circ_id, region.region_id, RBP_data.RBP_name, coverage, RBP_dist, sep='\t', file=out)

        else:
            print(region.circ_id, region.region_id, '.', coverage, '.', sep='\t', file=out)
