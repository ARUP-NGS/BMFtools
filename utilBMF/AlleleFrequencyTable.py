#!/usr/bin/env python

import argparse

from utilBMF.HTSUtils import printlog as pl
from MawCluster.PileupUtils import AlleleFrequenciesByBase


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
        "-p",
        "--progRepInterval",
        help="Number of positions between progress reports.",
        default=10000,
        type=int
        )
    parser.add_argument(
        "-b",
        "--bed",
        help="Path to bedfile.",
        default="default"
        )
    parser.add_argument(
        "--minMQ",
        "-m",
        help="Minimum mapping quality for inclusion. Default: 0.",
        default=0,
        type=int)
    parser.add_argument(
        "--minBQ",
        "-q",
        help="Minimum base quality score for inclusion. Default: 20.",
        default=20,
        type=int)
    args = parser.parse_args()
    minMQ = args.minMQ
    minBQ = args.minBQ
    if(args.bed == "default"):
        Output = AlleleFrequenciesByBase(
            args.bam,
            outputTsv=args.AlleleFreqTsv,
            progRepInterval=int(args.progRepInterval),
            minMQ=minMQ,
            minBQ=minBQ)
    else:
        Output = AlleleFrequenciesByBase(
            args.bam,
            outputTsv=args.AlleleFreqTsv,
            progRepInterval=int(args.progRepInterval),
            bedfile=args.bed,
            minMQ=minMQ,
            minBQ=minBQ)
    return Output


if(__name__ == "__main__"):
    main()
