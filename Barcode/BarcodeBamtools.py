#!/mounts/anaconda/bin/python

import pysam
from pysam import Samfile
from Bio import SeqIO
import argparse


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--fq', help="Provide your fastq file(s), post-processing with BarcodeUtils.py. If paired-end, provide read1 file first, followed by read2", nargs = "+", metavar = ('reads'),required=True)
    parser.add_argument('-p','--paired-end',help="Whether the experiment is paired-end or not. Default: True",default=True)
    parser.add_argument('-b','--bam-file',help="Aligned bam file, ready for the addition of tags.",required=True)
    parser.add_argument('-o','--output',help="Output path for finished BAM file. Defaults to basename.tagged.bam",default="default")
    args=parser.parse_args()
    if(args.paired_end==False):
        #I will write in how this is done later. Right now I am working on the more interesting and useful problem of paired-end
        print("I will code in this alternate version later. Right now I am working on the more interesting and useful problem of paired-end. Sorry!")
        return 
    if(len(args.fq) != 2):
        raise NameError("I DON'T KNOW WHAT WE'RE YELLING ABOUT. REQUIRED: TWO FASTQ FILES")
    read1 = SeqIO.parse(args.fq[0], "fastq")
    read2 = SeqIO.parse(args.fq[1], "fastq")
    inBAM = removeSecondary(args.bam_file)
    postFilterBAM = Samfile(inBAM,"rb")
    output=args.output
    if(output=="default"):
        output='.'.join(args.bam_file.split('.')[0:-1])+'.tagged.bam'
    outBAM = Samfile(output,"wb",template=postFilterBAM)
    for entry in postFilterBAM:
        if(entry.is_read1):
            #print("Read is read 1. Now getting variables from read 1's files for the tags.")
            tempRead = read1.next()
        elif(entry.is_read2):
            #print("Read is read 2. Now getting variables from read 2's files for the tags.")
            tempRead = read2.next()
        descArray=tempRead.description.split("###")
        entry.tags = entry.tags + [("BS",descArray[-2].strip())]
        entry.tags = entry.tags + [("FM",descArray[-1].strip())]
        if(descArray[-3].strip() == "AdapterPass"):
            entry.tags = entry.tags + [("AL",1)]
        else:
            entry.tags = entry.tags + [("AL",0)]
        outBAM.write(entry)
    outBAM.close()
    postFilterBAM.close()
    return

def removeSecondary(inBAM,outBAM="default"):
    from subprocess import call
    if(outBAM=="default"):
        outBAM = '.'.join(inBAM.split('.')[0:-1])+'.2ndrm.bam'
    command = "samtools view -h -F {} -o {}".format(inBAM,outBAM)
    call(command,shell=True)
    return outBAM

if(__name__=="__main__"):
    main()
