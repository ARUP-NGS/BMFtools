#!/usr/bin/env python
import argparse
import sys

from MawCluster.VCFWriters import SNVCrawler
from MawCluster.BCVCF import VCFStats
from MawCluster import BCVCF
from MawCluster import BCFastq
from BMFMain.ProcessingSteps import pairedFastqShades
from utilBMF import HTSUtils
from MawCluster.TLC import BMFXLC as CallIntraTrans
from pudb import set_trace

"""
bmftools contains various utilities for barcoded reads and for
somatic variant calling. Written to be in similar form to bcftools
and samtools.
"""


def main():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="bmfsuites")
    VCFStatsParser = subparsers.add_parser("vcfstats",
                                           description="Gets counts and"
                                           " frequencies for all SNV tr"
                                           "ansitions.")
    VCFHetParser = subparsers.add_parser("findhets",
                                         description="Split VCF lines with mul"
                                         "tiple ALTs, then writes to files onl"
                                         "y those variants with either desired"
                                         " Het frequency or given minor allele"
                                         " frequency.")
    DMultiPlexParser = subparsers.add_parser("dmp",
                                             description="Marks, combines, and"
                                             " processes a dataset of fastqs f"
                                             "or further analysis.")
    BamTagsParser = subparsers.add_parser("tagbam",
                                          description="Tags a BAM file with in"
                                          "formation from Fastq file(s)")
    SNVParser = subparsers.add_parser("snv", description="Call SNVs. Assumes "
                                      "that reads have been collapsed from a "
                                      "family size of at least 2")
    SVParser = subparsers.add_parser("sv",
                                     description="Call structural variants. R"
                                     "equires an Input BAM, coordinate-sorted"
                                     " and indexed, with BMF SV Tags included"
                                     ", and a BED File.")
    SMAParser = subparsers.add_parser(
        "sma",
        description="Tool for splitting a VCF File with multiple alts per"
        " line into a VCF where each line has a unique alt.")
    SNVParser.add_argument("inBAM",
                           help="Input BAM, Coordinate-sorted and indexed, "
                           "with BMF Tags included. bmftools runs on unflat"
                           "tened BAMs, but the results are not reliable be"
                           "cause this program assumes as much.")
    SNVParser.add_argument(
        "--bed",
        "-b",
        help="Full path to bed file.",
        default="default",
        metavar="bedpath")
    SNVParser.add_argument(
        "-o",
        "--outVCF",
        help="Output VCF File.",
        default="default",
        metavar="OutputVCF")
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
        default=10,
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
                          " for clustering families together",
                          default=35)
    SMAParser.add_argument(
        "inVCF",
        help="Input VCF", type=str)
    SMAParser.add_argument(
        "--outVCF",
        "-o",
        help="Output VCF. If unset, defaults to a modified form of the input.",
        default="default")
    VCFHetParser.add_argument(
        "inVCF", help="Input VCF", type=str)
    VCFHetParser.add_argument(
        "--outVCF",
        "-o",
        help="Output VCF. If unset, defaults to a modified form of the input.",
        default="default")
    VCFHetParser.add_argument(
        "--vcf-format",
        "-f",
        help="VCF Format. 'ExAC' for ExAC style, 'UK10K' for UK10K style.",
        default="ExAC")
    VCFHetParser.add_argument(
        "--min-het-frac",
        help="Minimum fraction of population het OR minimum minor "
        "allele frequency, depending on dataset.",
        type=float,
        default=0.025)
    # set_trace()

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
            bedFilteredVCF = BCVCF.FilterVCFFileByBed(OutVCF, bedfile=args.bed)
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
    if(args.bmfsuites == "sma"):
        Output = BCVCF.ISplitMultipleAlts(args.inVCF, outVCF=args.outVCF)
    if(args.bmfsuites == "findhets"):
        smaOut = BCVCF.ISplitMultipleAlts(args.inVCF, outVCF=args.outVCF)
        print("Multiple alts split: {}".format(smaOut))
        if(args.vcf_format.lower() == "exac"):
            hetOut = BCVCF.GetPotentialHetsVCF(smaOut,
                                               minHetFrac=args.min_het_frac,
                                               outVCF=args.outVCF)
        elif(args.vcf_format.lower() == "uk10k"):
            hetOut = BCVCF.GetPotentialHetsVCFUK10K(args.inVCF)
        print("Potential Hets VCF: {}".format(hetOut))
        return hetOut
    return 0


if(__name__ == "__main__"):
    main()
