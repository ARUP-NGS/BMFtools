#!/usr/bin/env python

import argparse

from BMFUtils.HTSUtils import printlog as pl
from MawCluster.BCVCF import AlleleFrequenciesByBase


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--bam',
        '-i',
        help=("Coordinate-Sorted, Indexed Bam File"),
        required=True
        )
    parser.add_argument(
        "-o",
        "--AlleleFreqTsv",
        help="Output Allele Frequency Tsv",
        default="default")
    parser.add_argument(
        "--p",
        "--progRepInterval",
        help="Number of positions between progress reports.",
        default=10000
        )
    args = parser.parse_args()
    Output = AlleleFrequenciesByBase(
        args.bam,
        outputTsv=args.AlleleFreqTsv,
        progRepInterval=int(args.progRepInterval))
    return Output


if(__name__ == "__main__"):
    main()
