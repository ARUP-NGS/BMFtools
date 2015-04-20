#!/usr/bin/env python
import argparse
import sys
import warnings
import cStringIO

from MawCluster.VCFWriters import SNVCrawler
from MawCluster.BCVCF import (VCFStats, CheckStdCallsForVCFCalls,
                              CheckVCFForStdCalls)
from MawCluster import BCVCF
from MawCluster import BCFastq
from BMFMain.ProcessingSteps import pairedFastqShades
from utilBMF import HTSUtils
from MawCluster.TLC import BMFXLC as CallIntraTrans
from MawCluster.FFPE import TrainAndFilter, FilterByDeaminationFreq
#  from pudb import set_trace

"""
bmftools contains various utilities for barcoded reads and for
somatic variant calling. Written to be similar in form to bcftools
and samtools.
"""

#  warnings.filterwarnings('error')


def main():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="bmfsuites")
    FFPEParser = subparsers.add_parser(
        "ffpe", description="Estimates cytosine deamination frequencies in a "
        "VCF outputs a VCF with appropriate variants filtered. FILTER: Deamin"
        "ationNoise")
    VCFCmpParser = subparsers.add_parser(
        "vcfcmp", description="Compares VCF files.")
    VCFStatsParser = subparsers.add_parser("vcfstats",
                                           description="Gets counts and"
                                           " frequencies for all SNV tr"
                                           "ansitions.")
    VCFAFFilterParser = subparsers.add_parser("faf",
                                              description="Filter a VCF by AF")
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
        default=80,
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
    SNVParser.add_argument("-m", "--min-frac-agreed",
                           help="Minimum fraction of family agreed on base.",
                           default=0.75, type=float)
    SNVParser.add_argument("--minFA", help="Minimum family agreed on base.",
                           default=3, type=int)
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
    VCFAFFilterParser.add_argument(
        "inVCF", help="Path to input VCF")
    VCFAFFilterParser.add_argument(
        "--maxAF", help="Max allele frequency to keep",
        type=float, default=0.1)
    VCFAFFilterParser.add_argument(
        "--outVCF", help="output VCF. 'default' will pick a modified form of"
        " the inVCF name.",
        type=str, default="default")
    FFPEParser.add_argument(
        "inVCF", help="Path to the input VCF.", type=str)
    FFPEParser.add_argument(
        "--pVal", help="P value for confidence interval. Default: 0.05",
        default=0.05, type=float)
    FFPEParser.add_argument(
        "--ctfreq", help="Estimated deamination frequency for the sample.",
        type=float, default=-1.)
    FFPEParser.add_argument(
        "--maxFreq", help="Maximum frequency for a C-[TU]/G-A event to be as"
        "sumed to be deamination due to formalin fixation.", default=0.1,
        type=float)
    VCFCmpParser.add_argument(
        "queryVCF", help="Query VCF to compare to reference VCF.",
        type=str)
    VCFCmpParser.add_argument(
        "--std", help="Reference VCF for comparing to query VCF.",
        type=str, required=True)
    VCFCmpParser.add_argument(
        "-o", "--outfile", help="Set output file path instead of stdout.",
        default=None)
    VCFCmpParser.add_argument(
        "--check-std",
        help=("If set, check standard VCF for calls in the query VCF rather "
              "than default behavior, which is checking the query VCF for "
              "calls that should be in the standard."),
        action="store_true", default=False)
    # set_trace()

    args = parser.parse_args()
    commandStr = " ".join(sys.argv)
    if(args.bmfsuites == "snv"):
        """
        import cProfile
        import pstats
        pr = cProfile.Profile()
        pr.enable()
        """
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
                                OutVCF=args.outVCF,
                                minFracAgreed=args.min_frac_agreed,
                                minFA=args.minFA)
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
        """
        s = cStringIO.StringIO()
        pr.disable()
        sortby = "cumulative"
        ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
        ps.print_stats()
        open("cProfile.stats.txt", "w").write(s.getvalue())
        """
        sys.exit(0)
    if(args.bmfsuites == "vcfstats"):
        OutTable = VCFStats(args.inVCF)
        sys.exit(0)
    if(args.bmfsuites == "dmp"):
        OutFastq1, OutFastq2 = pairedFastqShades(
            args.inFqs[0],
            args.inFqs[1],
            indexfq=args.indexFq)
        sys.exit(0)
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
        sys.exit(0)
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
    if(args.bmfsuites == "faf"):
        Output = BCVCF.IFilterByAF(args.inVCF, maxAF=args.maxAF,
                                   outVCF=args.outVCF)
        sys.exit(0)
    if(args.bmfsuites == "ffpe"):
        # ctfreq defaults to -1. This conditional checks to see if it's set.
        if(args.ctfreq < 0):
            Output = TrainAndFilter(args.inVCF, maxFreq=args.maxFreq,
                                    pVal=args.pVal)
        else:
            Output = FilterByDeaminationFreq(args.inVCF, pVal=args.pVal,
                                             ctfreq=args.ctfreq)
        sys.exit(0)
    if(args.bmfsuites == "vcfcmp"):
        if(args.check_std):
            if(args.outfile is not None):
                CheckStdCallsForVCFCalls(args.queryVCF, std=args.std,
                                         outfile=args.outfile)
            else:
                CheckStdCallsForVCFCalls(args.queryVCF, std=args.std,
                                         outfile=sys.stdout)
            sys.exit(0)
        if(args.outfile is not None):
            CheckVCFForStdCalls(args.queryVCF, std=args.std,
                                outfile=args.outfile)
        else:
            CheckVCFForStdCalls(args.queryVCF, std=args.std,
                                outfile=sys.stdout)
    sys.exit(0)


if(__name__ == "__main__"):
    main()
