import argparse
import logging

from Bio import SeqIO
import NewBarcodeSteps
# Contains utilities for the completion of a variety of
# tasks related to barcoded protocols for ultra-low
# frequency variant detection, particularly for circulating tumor DNA
#
# Sept. 15, 2014: Redesigned. New homing sequence: GACGG(T)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-i',
        '--fq',
        help="Provide your fastq file(s).",
        nargs="+",
        metavar=('reads'))
    parser.add_argument(
        '-r',
        '--ref',
        help="Prefix for the reference index. Required.",
        required=True)
    parser.add_argument(
        '--homing',
        help="Homing sequence for samples. If not set, defaults to CAGT.",
        metavar=('HomingSequence'),
        default="CAGT")
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
    if(args.logfile != "default"):
        logfile = args.logfile
    else:
        logfile = args.fq[0].split('.')[0] + '.log'
    logging.basicConfig(filename=logfile, level=logging.INFO)
    aligner, homing = args.aligner, args.homing
    ref, opts, bed = args.ref, args.opts, args.bed
    logging.info("Paired-end: {}".format(args.paired_end))
    if(args.paired_end is False or args.paired_end.lower() == "false"):
        if(args.initialStep == 1):
            logging.info("Beginning fastq processing.")
            consFq = NewBarcodeSteps.singleFastqProc(args.fq[0], homing=homing)
            logging.info("Beginning BAM processing.")
            TaggedBam = NewBarcodeSteps.singleBamProc(
                consFq,
                ref,
                opts,
                aligner=aligner)
            logging.info("Beginning VCF processing.")
            CleanParsedVCF = NewBarcodeSteps.singleVCFProc(TaggedBam, bed, ref)
            logging.info(
                "Last stop! Watch your step.")
            return
        elif(args.initialStep == 2):
            consFq = args.fq[0]
            logging.info("Beginning BAM processing.")
            TaggedBam = NewBarcodeSteps.singleBamProc(
                consFq,
                ref,
                opts,
                aligner=aligner)
            logging.info("Beginning VCF processing.")
            CleanParsedVCF = NewBarcodeSteps.singleVCFProc(TaggedBam, bed, ref)
            logging.info(
                "Last stop! Watch your step.")
            return
        elif(args.initialStep == 3):
            ConsensusBam = args.BAM
            logging.info("Beginning VCF processing.")
            CleanParsedVCF = NewBarcodeSteps.singleVCFProc(
                ConsensusBam,
                bed,
                ref)
            logging.info(
                "Last stop! Watch your step.")
            return
        else:
            raise ValueError("You have chosen an illegal initial step.")
    elif(args.paired_end or args.paired_end.lower() == "true"):
        if(args.initialStep == 1):
            logging.info("Beginning fastq processing.")
            trimfq1, trimfq2, trimfqSingle = NewBarcodeSteps.pairedFastqProc(
                args.fq[0], args.fq[1], homing=homing)
            logging.info("Beginning BAM processing.")
            procSortedBam = NewBarcodeSteps.pairedBamProc(
                trimfq1,
                trimfq2,
                consfqSingle=trimfqSingle,
                aligner=aligner,
                ref=ref)
            CleanParsedVCF = NewBarcodeSteps.pairedVCFProc(
                procSortedBam,
                ref=ref,
                opts=opts,
                bed=bed)
            logging.info(
                "Last stop! Watch your step.")
        elif(args.initialStep == 2):
            logging.info("Beginning BAM processing.")
            procSortedBam = NewBarcodeSteps.pairedBamProc(
                trimfq1,
                trimfq2,
                consfqSingle="default",
                aligner=aligner,
                ref=ref)
            logging.info("Beginning VCF processing.")
            CleanParsedVCF = NewBarcodeSteps.pairedVCFProc(
                procSortedBam,
                ref=ref,
                opts=opts,
                bed=bed)
            logging.info(
                "Last stop! Watch your step")
        elif(args.initialStep == 3):
            logging.info("Beginning VCF processing.")
            CleanParsedVCF = NewBarcodeSteps.pairedVCFProc(
                args.BAM,
                ref=ref,
                opts=opts,
                bed=bed)
            logging.info(
                "Last stop! Watch your step.")
        return

    else:
        raise ValueError(
            "Not a valid Boolean value for paired-end.")

if(__name__ == "__main__"):
    main()
