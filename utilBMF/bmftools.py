#!/usr/bin/env python
import argparse
import sys

from MawCluster.VCFWriters import SNVCrawler
"""
bmftools contains various utilities for barcoded reads and for
somatic variant calling. Written to be in similar form to bcftools
and samtools.
"""


def main():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="bmfsuites")
    SNVParser = subparsers.add_parser("snv")
    SNVParser.add_argument("inBAM",
                           help="Input BAM, Coordinate-sorted and indexed, "
                           "with BMF Tags included. bmftools runs on unflat"
                           "tened BAMs, but the results are not reliable be"
                           "cause this program assumes as much.")
    SNVParser.add_argument(
        "--bed",
        "-b",
        help="Full path to bed file.",
        default="default")
    SNVParser.add_argument(
        "--minBQ",
        help="Minimum Base Quality to consider",
        default=30,
        type=int)
    SNVParser.add_argument(
        "--minMQ",
        help="Minimum Mapping Quality to consider",
        default=10,
        type=int)
    SNVParser.add_argument(
        "--MaxPValue",
        "-p",
        help="Maximum P value to consider, in e notation.",
        type=float,
        default=1e-15)
    SNVParser.add_argument(
        "--keepConsensus",
        "-k",
        action="store_true")
    SNVParser.add_argument(
        "--logfile",
        help="Name for logfile.",
        default="default")
    SNVParser.add_argument(
        "-r",
        "--reference-fasta",
        help="Provide reference fasta.",
        required=True)
    args = parser.parse_args()
    commandStr = " ".join(sys.argv)

    if(args.bmfsuites == "snv"):
        # Beckon the Kraken
        if(args.bed != "default"):
            OutVCF = SNVCrawler(args.inBAM,
                                bed=args.bed,
                                minMQ=args.minMQ,
                                minBQ=args.minBQ,
                                MaxPValue=args.MaxPValue,
                                keepConsensus=args.keepConsensus,
                                commandStr=commandStr,
                                reference=args.reference_fasta,
                                reference_is_path=True)
        else:
            OutVCF = SNVCrawler(args.inBAM,
                                minMQ=args.minMQ,
                                minBQ=args.minBQ,
                                MaxPValue=args.MaxPValue,
                                keepConsensus=args.keepConsensus,
                                commandStr=commandStr,
                                reference=args.reference_fasta,
                                reference_is_path=True)
    return


if(__name__ == "__main__"):
    main()
