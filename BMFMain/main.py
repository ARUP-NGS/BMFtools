#!/usr/bin/env python
# cython: boundscheck=False
# cython: profile=True
# cython: cdivision=True
# cython: cdivision_warnings=True
import argparse
import logging
import os
import os.path
import pudb
import sys
import datetime
import subprocess

import BMFMain.ProcessingSteps as ps
from utilBMF.HTSUtils import printlog as pl, FacePalm, ThisIsMadness
from utilBMF import HTSUtils

"""
Contains utilities for the completion of a variety of
tasks related to barcoded protocols for ultra-low
frequency variant detection, particularly for circulating tumor DNA
Structural Variant detection tools are in active development.
"""

# Global Variables
global Logger
Logger = logging.getLogger("Primarylogger")


def main(argv=None):
    if argv is not None:
        sys.argv = argv
    # pudb.set_trace()
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'fq',
        help="Provide your fastq file(s).",
        nargs="+",
        metavar=('reads'))
    parser.add_argument(
        "--idxFastq",
        "-i",
        help="Path to index fastq",
        metavar="indexFastq")
    parser.add_argument(
        '--conf',
        help="Path to config file with settings.",
        metavar="ConfigPath",
        default="/yggdrasil/workspace/BMFTools/demo/config.txt")
    parser.add_argument(
        '-s',
        '--single-end',
        help="Whether the experiment is single-end or not. Default: False",
        default=False,
        type=bool)
    parser.add_argument(
        '--homing',
        help="Homing sequence for samples. If not set, defaults to GACGG.",
        metavar=('HomingSequence'),
        default="GACGG")
    parser.add_argument(
        '--shades',
        help="Use flag if using the shades barcode method.",
        default=True,
        action="store_true")
    parser.add_argument(
        '-a',
        '--aligner',
        help="Provide your aligner. Default: bwa",
        nargs='?',
        metavar='aligner',
        default='mem')
    parser.add_argument(
        '-o',
        '--opts',
        help="Additional aligner opts. E.g.: --opts '-L 0' ",
        nargs='?',
        default='')
    parser.add_argument(
        '-b',
        '--BAM',
        help="BAM file, if alignment has already run.",
        default="default")
    parser.add_argument(
        '--bed',
        help="full path to bed file used for variant-calling steps."
             "Can be indicated by the config file.",
        metavar="BEDFile",
        default="default")
    parser.add_argument(
        '--initialStep',
        help="1: Fastq. 2: Bam. 3. VCF",
        default=1,
        type=int)
    parser.add_argument(
        '-l',
        '--logfile',
        help="To change default logfile location.",
        default="default")
    parser.add_argument(
        '-p',
        '--file-prefix',
        help="Set non-default prefix.",
        default="default")
    parser.add_argument(
        '--minMQ',
        help="Minimum mapping quality for variant call inclusion. "
             "Can be indicated by the config file.",
        default=20)
    parser.add_argument(
        '--minBQ',
        help="Minimum base quality for variant call inclusion. "
             "Can be indicated by the config file.",
        default=90)
    parser.add_argument(
        "--minCov",
        help="Minimum coverage for including a position"
        " in the BamToCoverageBed",
        default=5
        )
    parser.add_argument(
        "--ref",
        "-r",
        default="default",
        help="Path to reference index. Can be indicated by the config file.")
    parser.add_argument(
        "--abrapath",
        default="default",
        help="Path to abra jar. Can be indicated by the config file.")
    parser.add_argument(
        "--barcodeIndex",
        help="If starting with the BAM step, provide the "
             "path to your barcode index as created.")
    parser.add_argument(
        "--lighter",
        help="Whether or not to use Lighter for error correction.",
        default="default")
    parser.add_argument(
        "--kmer",
        help="Kmer for error correction.",
        type=int)
    parser.add_argument(
        "--captureSize",
        help="Size of capture in base pairs. Required for Lighter.",
        type=int)
    parser.add_argument(
        "--alpha", help="Alpha parameter for Lighter.", type=float)
    parser.add_argument(
        "--p3Seq", help="3' primer sequence for cutadapt.", default="default")
    parser.add_argument(
        "--p5Seq", help="5' primer sequence for cutadapt.", default="default")
    parser.add_argument(
        "--review-dir", help="Prefix for review directory, where important re"
        "sults files will be moved at the end of analysis.", default="default")
    parser.add_argument("--minFA", help="Minimum family members agreed on bas"
                        "e for inclusion in variant call", default=2, type=int)
    args = parser.parse_args()
    confDict = HTSUtils.parseConfig(args.conf)
    if("minMQ" in confDict.keys()):
        minMQ = int(confDict['minMQ'])
    else:
        minMQ = args.minMQ
    if("minBQ" in confDict.keys()):
        minBQ = int(confDict['minBQ'])
    else:
        minBQ = args.minBQ
    if("abrapath" in confDict.keys()):
        abrapath = confDict['abrapath']
    else:
        abrapath = args.abrapath
    if("minFA" in confDict.keys()):
        minFA = int(confDict['minFA'])
    else:
        minFA = args.minFA
    if("minFracAgreed" in confDict.keys()):
        minFracAgreed = float(confDict['minFracAgreed'])
    else:
        minFracAgreed = args.minFracAgreed
    global gatkpath
    if("gatkpath" in confDict.keys()):
        gatkpath = confDict['gatkpath']
    captureSize = None
    kmer = None
    alpha = None
    if("lighter" in confDict.keys()):
        lighter = (confDict['lighter'].lower() == "true")
    if "lighter" in args and args.lighter != "default":
        lighter = (args.lighter.lower() == "true")
    if("lighter" in locals()):
        if("kmer" in confDict.keys()):
            kmer = int(confDict["kmer"])
        if(args.kmer is not None):
            kmer = args.kmer
        if("captureSize" in confDict.keys()):
            captureSize = confDict['captureSize']
        if(args.captureSize is not None):
            captureSize = int(args.captureSize)
        if("alpha" in confDict.keys()):
            alpha = confDict['alpha']
        if(args.alpha is not None):
            alpha = args.alpha
        pl("captureSize: {}".format(captureSize))
        if(captureSize is None):
            raise ThisIsMadness("required captureSize variable not set!")
    dateStr = datetime.datetime.now().strftime("%Y-%b-%d,%H-%m-%S")
    global reviewdir
    reviewdir = ""
    makeReviewDir = (args.review_dir != "default")
    if(makeReviewDir is True):
        reviewdir = ".".join([args.review_dir, dateStr, "reviewdir"])
        if(os.path.isdir(reviewdir)):
            raise ThisIsMadness("Review directory exists - even with "
                                "the time stamp. Abort!")
        if(os.path.isfile(reviewdir)):
            raise ThisIsMadness("Not only is the review directory name"
                                " taken, but it's not even a folder?")
        os.mkdir(reviewdir)
    pl("Review directory: {}")
    # Begin logging
    if(args.logfile != "default"):
        logfile = args.logfile
    else:
        logfile = (os.getcwd() + "/" +
                   os.path.basename(args.fq[0]).split('.')[0] +
                   '.log')
    if(args.file_prefix != "default"):
        logfile = os.getcwd() + "/" + args.file_prefix + ".log"
    if(os.path.isfile(logfile)):
        os.remove(logfile)
        pl("Log file existed - deleting!")

    # Logger which holds both console and file loggers
    Logger.setLevel(logging.DEBUG)

    # Console handler - outputs to console.
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)

    # create formatter
    formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s')

    # add formatter to ch
    ch.setFormatter(formatter)

    # add ch to logger
    Logger.addHandler(ch)

    # File logger - outputs to log file.
    fl = logging.FileHandler(filename=logfile)
    fl.setFormatter(formatter)
    fl.setLevel(logging.DEBUG)

    Logger.addHandler(fl)

    pl("Log file is {}".format(logfile))
    pl("Command string to call BMFTools: python {}".format(" ".join(sys.argv)))
    pl("Note: You may need to add quotes for more complicated options.")
    if("aligner" in confDict.keys()):
        aligner = confDict['aligner']
    else:
        aligner = args.aligner
    homing = args.homing
    if("ref" in confDict.keys()):
        ref = confDict['ref']
    else:
        if(args.ref != "default"):
            ref = args.ref
        else:
            HTSUtils.FacePalm("Reference fasta required either in "
                              "config file or by command-line options.!")
    opts = args.opts
    if("opts" in confDict.keys()):
        opts = args.opts
    if("bed" in confDict.keys()):
        bed = confDict['bed']
    else:
        if(args.bed != "default"):
            bed = args.bed
        else:
            FacePalm("Bed file required for analysis.")
    if("single_end" in confDict.keys()):
        if(confDict['single_end'].lower() == "true"):
            single_end = True
        else:
            single_end = False
    else:
        single_end = args.single_end
    if("p3Seq" in confDict.keys()):
        p3Seq = confDict['p3Seq']
    if(args.p3Seq != "default"):
        p3Seq = args.p3Seq
    if("p5Seq" in confDict.keys()):
        p5Seq = confDict['p5Seq']
    if(args.p5Seq != "default"):
        p5Seq = args.p5Seq
    if("bwapath" in confDict.keys()):
        bwapath = confDict["bwapath"]
    if(single_end is True):
        pl("Single-end analysis chosen.")
        HTSUtils.ThisIsMadness("Single-end analysis not currently "
                               "supported. Soon!")
    else:
        pl("Paired-end analysis chosen.")
        if(args.initialStep == 1):
            pl("Beginning fastq processing.")
            if kmer is None and alpha is not None:
                trimfq1, trimfq2 = ps.pairedFastqShades(
                    args.fq[0], args.fq[1], indexfq=args.idxFastq,
                    lighter=lighter, captureSize=captureSize, alpha=alpha,
                    p3Seq=p3Seq, p5Seq=p5Seq)
            else:
                trimfq1, trimfq2 = ps.pairedFastqShades(
                    args.fq[0], args.fq[1], indexfq=args.idxFastq,
                    lighter=lighter, kmer=kmer, captureSize=captureSize,
                    p3Seq=p3Seq, p5Seq=p5Seq)
            if(makeReviewDir):
                subprocess.check_call(["cp", trimfq1, trimfq2, reviewdir])
            if("bwapath" in locals()):
                procSortedBam = ps.pairedBamProc(
                    trimfq1,
                    trimfq2,
                    aligner=aligner,
                    ref=ref,
                    bed=bed,
                    mincov=int(args.minCov),
                    abrapath=abrapath,
                    bwapath=bwapath)
            else:
                procSortedBam = ps.pairedBamProc(
                    trimfq1, trimfq2,
                    aligner=aligner, ref=ref,
                    bed=bed,
                    mincov=int(args.minCov), abrapath=abrapath)
            if(makeReviewDir):
                subprocess.check_call(["cp", procSortedBam, reviewdir])
            VCFOutDict = ps.pairedVCFProc(
                procSortedBam,
                ref=ref,
                opts=opts,
                bed=bed,
                minMQ=minMQ,
                minBQ=minBQ,
                commandStr=" ".join(sys.argv),
                minFA=minFA, minFracAgreed=minFracAgreed)
            for key in VCFOutDict.keys():
                if(makeReviewDir):
                    subprocess.check_call(["cp", VCFOutDict[key], reviewdir])
            pl("Last stop! Watch your step.")
        elif(args.initialStep == 2):
            pl("Beginning BAM processing.")
            if("bwapath" not in locals()):
                procSortedBam = ps.pairedBamProc(
                    args.fq[0],
                    args.fq[1],
                    aligner=aligner,
                    ref=ref,
                    bed=bed,
                    mincov=int(args.minCov),
                    abrapath=abrapath,
                    bwapath=bwapath, barIndex=args.barcodeIndex)
            else:
                procSortedBam = ps.pairedBamProc(
                    args.fq[0],
                    args.fq[1],
                    aligner=aligner,
                    ref=ref,
                    bed=bed,
                    mincov=int(args.minCov),
                    abrapath=abrapath,
                    barIndex=args.barcodeIndex)
            pl("Beginning VCF processing.")
            CleanParsedVCF = ps.pairedVCFProc(
                procSortedBam,
                ref=ref,
                opts=opts,
                bed=bed,
                minMQ=minMQ,
                minBQ=minBQ,
                commandStr=" ".join(sys.argv),
                minFA=minFA, minFracAgreed=minFracAgreed)
            pl("Last stop! Watch your step")
        elif(args.initialStep == 3):
            pl("Beginning VCF processing.")
            CleanParsedVCF = ps.pairedVCFProc(
                args.BAM,
                ref=ref,
                opts=opts,
                bed=bed,
                minMQ=args.minMQ,
                minBQ=args.minBQ,
                reference=args.ref,
                commandStr=" ".join(sys.argv),
                minFA=minFA, minFracAgreed=minFracAgreed)
            pl("Last stop! Watch your step.")
    if(makeReviewDir):
        subprocess.check_call(["cp", logfile, reviewdir])
    return

global __version__

__version__ = "0.0.7.0"

if(__name__ == "__main__"):
    sys.exit(main())
