#!/usr/bin/env python
import argparse
import sys

from MawCluster.VCFWriters import SNVCrawler
from MawCluster.BCVCF import VCFStats
from MawCluster import BCFastq
from BMFMain.ProcessingSteps import pairedFastqShades
from utilBMF import HTSUtils
from MawCluster.TLC import BMFXLC as CallIntraTrans

"""
bmftools contains various utilities for barcoded reads and for
somatic variant calling. Written to be in similar form to bcftools
and samtools.
"""


def main():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="bmfsuites")
    VCFStatsParser = subparsers.add_parser("vcfstats")
    DMultiPlexParser = subparsers.add_parser("dmp")
    BamTagsParser = subparsers.add_parser("tagbam")
    SNVParser = subparsers.add_parser("snv")
    SVParser = subparsers.add_parser("sv")
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
        "-o",
        "--outVCF",
        help="Output VCF File.",
        default="default")
    SNVParser.add_argument(
        "--minBQ",
        help="Minimum Base Quality to consider",
        default=30,
        type=int)
    SNVParser.add_argument(
        "--minMQ",
        help="Minimum Mapping Quality to consider",
        default=20,
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
    VCFStatsParser.add_argument(
        "inVCF",
        help="Input VCF, as created by SNVCrawler.")
    DMultiPlexParser.add_argument(
        "inFqs",
        nargs="+",
        help="Input Fastq Files")
    DMultiPlexParser.add_argument(
        "-i",
        "--indexFq",
        metavar="indexFastq",
        help="Index Fastq")
    BamTagsParser.add_argument(
        "inBAM",
        metavar="inBAM",
        help="Untagged Bam")
    BamTagsParser.add_argument(
        "--fastq",
        "-f",
        metavar="InFastqs",
        nargs="+",
        help="Tagged, Merged Fastq File")
    SVParser.add_argument(
        'bam',
        help=("Coordinate-Sorted, Indexed Bam File"),
        )
    SVParser.add_argument(
        "-b",
        "--bed",
        help="Path to bedfile.",
        default="default"
        )
    SVParser.add_argument(
        "--minMQ",
        "-m",
        help="Minimum mapping quality for inclusion. Default: 0.",
        default=0,
        type=int)
    SVParser.add_argument(
        "--minBQ",
        help="Minimum base quality for inclusion. Default: 0.",
        default=0,
        type=int)
    SVParser.add_argument(
        "-o",
        "--outTsv",
        help="Output tsv",
        default="default"
        )
    SVParser.add_argument(
        "--minPileupLen",
        "-l",
        help="Length of interval to be considered for call.",
        default=8,
        type=int)
    SVParser.add_argument(
        "--minClustDepth",
        "-d",
        default=10,
        type=int,
        help="Minimum depth for a cluster to be considered for call.")
    SVParser.add_argument(
        "--ref",
        "-r",
        help="Path to reference index.",
        required=True)
    SVParser.add_argument("--insert-distance",
                          "-i",
                          help="Maximum difference between edit distances"
                          " for clustering families together")

    args = parser.parse_args()
    commandStr = " ".join(sys.argv)

    if(args.bmfsuites == "snv"):
        # Beckon the Kraken
        if(args.bed != "default"):
            OutVCF = SNVCrawler(args.inBAM,
                                bed=args.bed,
                                minMQ=args.minMQ,
                                minBQ=args.minBQ,
                                MaxPValue=float(args.MaxPValue),
                                keepConsensus=args.keepConsensus,
                                commandStr=commandStr,
                                reference=args.reference_fasta,
                                reference_is_path=True,
                                OutVCF=args.outVCF)
            OutTable = VCFStats(OutVCF)
        else:
            OutVCF = SNVCrawler(args.inBAM,
                                minMQ=args.minMQ,
                                minBQ=args.minBQ,
                                MaxPValue=args.MaxPValue,
                                keepConsensus=args.keepConsensus,
                                commandStr=commandStr,
                                reference=args.reference_fasta,
                                reference_is_path=True,
                                OutVCF=args.outVCF)
            OutTable = VCFStats(OutVCF)
    if(args.bmfsuites == "vcfstats"):
        OutTable = VCFStats(args.inVCF)
    if(args.bmfsuites == "dmp"):
        OutFastq1, OutFastq2 = pairedFastqShades(
            args.inFqs[0],
            args.inFqs[1],
            indexfq=args.indexFq)
    if(args.bmfsuites == "sv"):
        if(args.bed == "default"):
            HTSUtils.FacePalm("Bed file required!")
        else:
            Output = CallIntraTrans(
                args.bam,
                outfile=args.outTsv,
                bedfile=args.bed,
                minMQ=args.minMQ,
                minBQ=args.minBQ,
                minClustDepth=args.minClustDepth,
                minPileupLen=args.minPileupLen,
                ref=args.ref,
                insDistance=args.insert_distance)
        return Output
    return 0


if(__name__ == "__main__"):
    main()
