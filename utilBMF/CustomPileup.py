#!/usr/bin/env python

import argparse

from MawCluster.PileupUtils import CustomPileupToTsv


def main():
    parser = argparse.ArgumentParser()
    '''
                          PileupTsv="default",
                      TransitionTable="default",
                      StrandedTTable="default",
                      bedfile="default",
                      progRepInterval=1000
    '''
    parser.add_argument(
        '--bam',
        '-i',
        help=("Coordinate-Sorted, Indexed Bam File"),
        required=True
        )
    parser.add_argument(
        '--bed',
        '-b',
        help="Bedfile to use.",
        required=True)
    parser.add_argument(
        "-o",
        "--PileupTsv",
        help="Output Pileup Tsv",
        default="default")
    parser.add_argument(
        "-t",
        "--TransitionTable",
        help="Output transition table.",
        default="default")
    parser.add_argument(
        "-s",
        "--StrandedTTable",
        help="Output stranded transition table.",
        default="default")
    parser.add_argument(
        "-p",
        "--progRepInterval",
        help="Number of positions between progress reports.",
        default=10000
        )
    parser.add_argument(
        "--minMQ",
        help="Minimum mapping quality.",
        default=10
        )
    parser.add_argument(
        "--minBQ",
        help="Minimum base quality.",
        default=20
        )
    args = parser.parse_args()
    Output = CustomPileupToTsv(args.bam,
                               PileupTsv=args.PileupTsv,
                               TransitionTable=args.TransitionTable,
                               StrandedTTable=args.StrandedTTable,
                               bedfile=args.bed,
                               progRepInterval=int(args.progRepInterval),
                               minMQ=args.minMQ,
                               minBQ=args.minBQ)
    return Output


if(__name__ == "__main__"):
    main()
