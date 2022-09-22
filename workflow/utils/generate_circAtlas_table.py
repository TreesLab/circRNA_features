#! /usr/bin/env python

import argparse
import re


TITLES = (
    'event_id',
    'circAtlas_id',
    'host_gene_id',
    'host_gene',
    'circRNA type',
    'Multiple Conservation Score(MCS)',
    'MCS-detail(species)',
    'MCS-detail(tissue)',
    'MCS-detail(sample)',
    'Tissue Specifity Index',
    'Length (nt)',
    'Type length',
    'Algorithm',
    'circRNA type(intron)',
    'circRNA type(exon)',
    'circRNA type(intergenic)',
    'circRNA type(Non-Repeat)',
    'circRNA type(antisense)',
    'circRNA type(exon/intergenic)',
    'circRNA type(3-utr)',
    'circRNA type(5-utr)',
    'Type length(mature)',
    'Type length(Estimated)',
    'Type length(Unknown)',
    'Algorithm(CIRI2)',
    'Algorithm(DCC)',
    'Algorithm(find_circ)',
    'Algorithm(CIRCexplorer2)'
)


class CircAtlasTableParser:
    _circRNA_types = (
        'intron',
        'exon',
        'intergenic',
        'Non-Repeat',
        'antisense',
        'exon/intergenic',
        '3-utr',
        '5-utr'
    )
    _type_length = ('mature', 'Estimated', 'Unknown')
    _algorithms = ('CIRI2', 'DCC', 'find_circ', 'CIRCexplorer2')

    @classmethod
    def parse(cls, circAtlas_data):
        for line in circAtlas_data:
            data = line.rstrip('\n').split('\t')

            event_id = cls._get_event_id(data[2], data[3])
            circAtlas_id = data[1]
            host_gene_id = cls._get_host_gene_id(data[5])
            host_gene = cls._get_host_gene(data[1])
            circRNA_types_array = cls._get_circRNA_types_array(data[4])
            type_length_array = cls._get_type_length_array(data[12])
            algorithms_array = cls._get_algorithms_array(data[13])

            result = [
                event_id,
                circAtlas_id,
                host_gene_id,
                host_gene,
                data[4],
                *data[6:14],
                *circRNA_types_array,
                *type_length_array,
                *algorithms_array
            ]

            yield result

    @staticmethod
    def _get_event_id(position, strand):
        chr_, pos1, pos2 = re.split(r'[:|]', position)

        chr_ = re.sub(r'^chr', '', chr_)
        pos1, pos2 = sorted(map(int, [pos1, pos2]))

        if strand == "+":
            donor, acceptor = pos2, pos1
        elif strand == "-":
            donor, acceptor = pos1, pos2

        return f"{chr_}:{donor}|{acceptor}({strand})"

    @staticmethod
    def _get_host_gene_id(gene_ids):
        if gene_ids == 'n/a':
            return ''
        else:
            gene_ids_without_ver = re.sub(r'\.[0-9]+', '', gene_ids)
            return gene_ids_without_ver

    @staticmethod
    def _get_host_gene(circAtlas_id):
        m = re.search(r'[^-]+-(.+)_[0-9]+', circAtlas_id)
        if m:
            host_gene = m[1]
            if host_gene == 'intergenic':
                return ''
            else:
                return host_gene

    @classmethod
    def _get_circRNA_types_array(cls, circRNA_type):
        convertor = BinaryArrayConvertor(cls._circRNA_types)
        array = convertor.convert([circRNA_type])
        return array

    @classmethod
    def _get_type_length_array(cls, type_length):
        convertor = BinaryArrayConvertor(cls._type_length)
        array = convertor.convert([type_length])
        return array

    @classmethod
    def _get_algorithms_array(cls, algorithms):
        convertor = BinaryArrayConvertor(cls._algorithms)
        array = convertor.convert(algorithms.split(','))
        return array


class BinaryArrayConvertor:
    def __init__(self, columns):
        self.columns = columns

    def convert(self, list_):
        array = ['0'] * (len(self.columns))

        for item in list_:
            if item in self.columns:
                array[self.columns.index(item)] = '1'

        return array


def create_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('circAtlas_table_file', type=argparse.FileType('r'))

    return parser


if __name__ == "__main__":
    parser = create_parser()
    args = parser.parse_args()

    print(*TITLES, sep='\t')

    _ = args.circAtlas_table_file.readline()

    for result in CircAtlasTableParser.parse(args.circAtlas_table_file):
        print(*result, sep='\t')
