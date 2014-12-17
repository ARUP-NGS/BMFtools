import argparse
import logging
import os.path

import ProcessingSteps as ps
from HTSUtils import printlog as pl
import HTSUtils
# Contains utilities for the completion of a variety of
# tasks related to barcoded protocols for ultra-low
# frequency variant detection, particularly for circulating tumor DNA


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-i',
        '--fq',
        help="Provide your fastq file(s).",
        nargs="+",
        metavar=('reads'))
    parser.add_argument(
        '--conf',
        help="Path to config file with settings.",
        metavar="ConfigPath",
        default="default")
    parser.add_argument(
        '-r',
        '--ref',
        help="Prefix for the reference index. Required.",
        default="default")
    parser.add_argument(
        '--homing',
        help="Homing sequence for samples. If not set, defaults to GACGG.",
        metavar=('HomingSequence'),
        default="GACGG")
    parser.add_argument(
        '--shades',
        help="Set to true if using the shades barcode method.",
        default="false")
    parser.add_argument(
        '-p',
        '--paired-end',
        help="Whether the experiment is paired-end or not. Default: True",
        default="True")
    parser.add_argument(
        '-a',
        '--aligner',
        help="Provide your aligner. Default: bwa",
        nargs='?',
        metavar='aligner',
        default='bwa')
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
        help="full path to bed file used for variant-calling steps.",)
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
    args = parser.parse_args()
    # Begin logging
    if(args.logfile != "default"):
        logfile = args.logfile
    else:
        logfile = args.fq[0][0:-6].split('/')[-1] + '.log'
    if(os.path.isfile(logfile)):
        os.remove(logfile)
        pl("Log file existed - deleting!")
    pl("Log file is {}".format(logfile))
    logging.basicConfig(filename=logfile,
                        level=logging.INFO,
                        format="%(levelname)s [%(asctime)s]: %(message)s")
    import sys
    pl("Command string to call BMFTools: python {}".format(' '.join(sys.argv)))
    aligner, homing = args.aligner, args.homing
    if(args.ref != "default"):
        ref = args.ref
    opts, bed = args.opts, args.bed
    confDict = {}
    if(args.conf != "default"):
        confLines = [line.strip().split("=") for line in open(
                     args.conf, "r").readlines() if line[0] != "#"]
        for line in confLines:
            confDict[line[0]] = line[1]
        try:
            ref = confDict["ref"]
        except KeyError:
            ref = args.ref
            pl("Config file does not contain a ref key.")
            if(ref == "default"):
                HTSUtils.FacePalm("Reference required either "
                                  "as command line option or config file.")
    pl("Paired-end: {}".format(args.paired_end))
    if(args.paired_end is False or args.paired_end.lower() == "false"):
        if(args.initialStep == 1):
            pl("Beginning fastq processing.")
            consFq = ps.singleFastqProc(args.fq[0], homing=homing)
            pl("Beginning BAM processing.")
            TaggedBam = ps.singleBamProc(
                consFq,
                ref,
                opts,
                aligner=aligner)
            pl("Beginning VCF processing.")
            CleanParsedVCF = ps.singleVCFProc(TaggedBam, bed, ref)
            pl("Last stop! Watch your step.")
            return
        elif(args.initialStep == 2):
            consFq = args.fq[0]
            pl("Beginning BAM processing.")
            TaggedBam = ps.singleBamProc(
                consFq,
                ref,
                opts,
                aligner=aligner)
            pl("Beginning VCF processing.")
            CleanParsedVCF = ps.singleVCFProc(TaggedBam, bed, ref)
            pl("Last stop! Watch your step.")
            return
        elif(args.initialStep == 3):
            ConsensusBam = args.BAM
            pl("Beginning VCF processing.")
            CleanParsedVCF = ps.singleVCFProc(
                ConsensusBam,
                bed,
                ref)
            pl("Last stop! Watch your step.")
            return
        else:
            raise ValueError("You have chosen an illegal initial step.")
    elif(args.paired_end.lower() == "true"):
        if(args.initialStep == 1):
            pl("Beginning fastq processing.")
            if(args.shades.lower() != "true"):
                pl("About to run"
                    " ps.pairedFastqProc("
                    "{}, {}, homing={})".format(
                        args.fq[0], args.fq[1], homing))
                (trimfq1, trimfq2, trimfqSingle,
                 barcodeIndex) = ps.pairedFastqProc(
                    args.fq[0], args.fq[1], homing=homing)
                pl("Beginning BAM processing.")
                procSortedBam = ps.pairedBamProc(
                    trimfq1,
                    trimfq2,
                    consfqSingle=trimfqSingle,
                    aligner=aligner,
                    ref=ref,
                    barIndex=barcodeIndex,
                    bed=bed)
            elif(args.shades.lower() == "true"):
                trimfq1, trimfq2, barcodeIndex = ps.pairedFastqShades(
                    args.fq[0], args.fq[1], args.fq[2])
                procSortedBam = ps.pairedBamProc(
                    trimfq1,
                    trimfq2,
                    aligner=aligner,
                    ref=ref,
                    barIndex=barcodeIndex,
                    bed=bed)
            CleanParsedVCF = ps.pairedVCFProc(
                procSortedBam,
                ref=ref,
                opts=opts,
                bed=bed)
            pl(
                "Last stop! Watch your step.")
        elif(args.initialStep == 2):
            pl("Beginning BAM processing.")
            procSortedBam = ps.pairedBamProc(
                trimfq1,
                trimfq2,
                consfqSingle="default",
                aligner=aligner,
                ref=ref)
            pl("Beginning VCF processing.")
            CleanParsedVCF = ps.pairedVCFProc(
                procSortedBam,
                ref=ref,
                opts=opts,
                bed=bed)
            pl(
                "Last stop! Watch your step")
        elif(args.initialStep == 3):
            pl("Beginning VCF processing.")
            CleanParsedVCF = ps.pairedVCFProc(
                args.BAM,
                ref=ref,
                opts=opts,
                bed=bed)
            pl("Last stop! Watch your step.")
        return

    else:
        raise ValueError(
            "Not a valid Boolean value for paired-end.")

if(__name__ == "__main__"):
    main()
