#! /usr/bin/env python

import sys
from itertools import groupby
from sam import SamFormat


def get_uniq_matches():
    sam_data_reader = map(SamFormat, sys.stdin)

    for qname, sam_gp in groupby(sam_data_reader, key=lambda data: data.qname):

        if qname is None:
            for sam_data in sam_gp:
                print(sam_data)
        else:
            sam_gp = list(sam_gp)
            flags = (sam_data.flag for sam_data in sam_gp)
            for flag in flags:
                if flag & 2048:
                    break
            else:
                for sam_data in sam_gp:
                    print(sam_data)


if __name__ == "__main__":
    get_uniq_matches()
