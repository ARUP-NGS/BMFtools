#!/usr/bin/env python

import argparse

from utilBMF.HTSUtils import printlog as pl, FacePalm
from MawCluster.TLC import XLocIntrachromosomalFusionCaller as CallTrans


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'bam',
        help=("Coordinate-Sorted, Indexed Bam File"),
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
        help="Minimum base quality for inclusion. Default: 0.",
        default=0,
        type=int)
    parser.add_argument(
        "-o",
        "--outTsv",
        help="Output tsv",
        default="default"
        )
    parser.add_argument(
        "--minPileupLen",
        "-l",
        help="Length of interval to be considered for call.",
        default=10,
        type=int)
    parser.add_argument(
        "--minClustDepth",
        "-d",
        default=5,
        type=int,
        help="Minimum depth for a cluster to be considered for call.")
    parser.add_argument(
        "--ref",
        "-r",
        help="Path to reference index.",
        required=True)
    parser.add_argument("--insert-distance",
                        "-i",
                        help="Maximum difference between edit distances"
                        " for clustering families together")
    args = parser.parse_args()
    if(args.bed == "default"):
        FacePalm("Bed file required!")
    else:
        Output = CallTrans(
            args.bam,
            outfile=args.outTsv,
            bedfile=args.bed,
            minMQ=args.minMQ,
            minBQ=args.minBQ,
            minClustDepth=args.minClustDepth,
            minPileupLen=args.minPileupLen,
            ref=args.ref)
    return Output


if(__name__ == "__main__"):
    main()
