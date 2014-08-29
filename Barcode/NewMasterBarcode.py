#/mounts/anaconda/bin/python

from Bio import SeqIO
import argparse
import pysam
import subprocess
import NewBarcodeSteps
#Contains utilities for the completion of a variety of 
#tasks related to barcoded protocols for ultra-low
#frequency variant detection, particularly for circulating tumor DNA
#
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--fq', help="Provide your fastq file(s). Note: if '--BAM' option is set, this needs to be the Family Fastq files", nargs = "+", metavar = ('reads'))
    parser.add_argument('-r','--ref',help="Prefix for the reference index. Required.",required=True)
    parser.add_argument('--homing', help="Homing sequence for samples. If not set, defaults to CAGT.",metavar=('HomingSequence'),default="CAGT")
    parser.add_argument('-p','--paired-end',help="Whether the experiment is paired-end or not. Default: True",default="True")
    parser.add_argument('-a','--aligner', help="Provide your aligner. E.g., 'bwa', 'bowtie2', or 'snap'. Currently only 'bwa' is supported. Default: bwa", nargs='?', metavar='aligner', default='bwa')
    parser.add_argument('-o','--opts', help="Please place additional arguments for the aligner after this tag in quotation marks. E.g.: --opts '-L 0' ", nargs='?', default='')
    parser.add_argument('-b','--BAM', help="BAM file, if alignment has already run and fastq family files created, or if the consensus BAM has been created.",default="default")
    parser.add_argument('--bed',help="full path to bed file used for variant-calling steps.",default="/yggdrasil/workspace/Barcode_Noah/WorkDir/xgen-pan-cancer-sorted.trim.bed")
    parser.add_argument('--initialStep',help="1: start with fastq's. 2: start with aligned bam file and the tagged fastq. 3: start with the final processes BAM for variant calling.",default=1,type=int) 
    args=parser.parse_args()
    aligner, homing=args.aligner, args.homing
    ref, opts, bed = args.ref, args.opts, args.bed
    print(args.paired_end)
    if(args.paired_end==False or args.paired_end.lower()=="false"):
        if(args.initialStep==1):
            consFq = NewBarcodeSteps.singleFastqProc(args.fq[0], ref,aligner=aligner, homing=homing,opts=opts)
            TaggedBam = NewBarcodeSteps.singleBamProc(consFq, ref, opts, aligner=aligner)
            CleanParsedVCF = NewBarcodeSteps.singleVCFProc(TaggedBam, bed, ref)
            print("This is far as the program goes at this point. Thank you for playing!")
            return
        elif(args.initialStep==2):
            outbam=args.BAM
            consFq = args.fq[0]
            TaggedBam = NewBarcodeSteps.singleBamProc(consFq, ref, opts, aligner=aligner)
            CleanParsedVCF = NewBarcodeSteps.singleVCFProc(TaggedBam, bed, ref)
            print("This is far as the program goes at this point. Thank you for playing!")
            return
        elif(args.initialStep==3):
            ConsensusBam = args.BAM
            CleanParsedVCF = NewBarcodeSteps.singleVCFProc(TaggedBam, bed, ref)
            print("This is far as the program goes at this point. Thank you for playing!")
            return
        else:
            raise ValueError("You have chosen an illegal initial step.")
    elif(args.paired_end==True or args.paired_end.lower()=="true"):
        if(args.initialStep==1):
            outbam, trimfq1, trimfq2 = NewBarcodeSteps.pairedFastqProcAlign(args.fq[0], args.fq[1], homing, ref=ref, aligner=aligner)
            procSortedBam = NewBarcodeSteps.pairedBamProc(trimfq1, trimfq2, outbam)
            CleanParsedVCF = NewBarcodeSteps.pairedVCFProc(procSortedBam, ref=ref, opts=opts, bed=bed)
            print("This is far as the program goes at this point. Thank you for playing!")
        elif(args.initialStep==2):
            procSortedBam = NewBarcodeSteps.pairedBamProc(args.fq[0], args.fq[1], args.BAM)
            CleanParsedVCF = NewBarcodeSteps.pairedVCFProc(procSortedBam, ref=ref, opts=opts, bed=bed)
            print("This is far as the program goes at this point. Thank you for playing!")
        elif(args.initialStep==3):
            CleanParsedVCF = NewBarcodeSteps.pairedVCFProc(args.BAM, ref=ref,opts=opts,bed=bed)    
            print("This is far as the program goes at this point. Thank you for playing!")
        return

    else:
        raise ValueError("You have not selected an appropriate Boolean value for paired-end.") 

if(__name__=="__main__"):
    main()
