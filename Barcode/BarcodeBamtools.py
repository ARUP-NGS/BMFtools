#!/mounts/anaconda/bin/python

import pysam
from pysam import Samfile
from Bio import SeqIO
from BarcodeUtils import main as part1
import argparse


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--fq', help="Provide your fastq file(s), post-processing with BarcodeUtils.py. If paired-end, provide read1 file first, followed by read2", nargs = "+", metavar = ('reads'),required=True)
    parser.add_argument('-p','--paired-end',help="Whether the experiment is paired-end or not. Default: True",default=True)
    parser.add_argument('-b','--bam-file',help="Aligned bam file, ready for the addition of tags.",required=True)
    parser.add_argument('-o','--output',help="Output path for finished BAM file. Defaults to basename.tagged.bam",default="default")
    args=parser.parse_args()
    if(args.paired_end==False or args.paired_end=="False"):
        reads = SeqIO.parse(args.fq[0],"fastq")
        inBAM = Samfile(args.bam_file,"rb")
        output=args.output
        if(output=="default"):
            output='.'.join(args.bam_file.split('.')[0:-1])+'.tagged.bam'
        outBAM = Samfile(output,"wb",template=inBAM)
        for entry in inBAM:
            tempRead = reads.next()
            descArray=tempRead.description.split("###")
            entry.tags = entry.tags + [("BS",descArray[-2].strip())]
            entry.tags = entry.tags + [("FM",descArray[-1].strip())]
            if(descArray[-3].strip() == "AdapterPass"):
                entry.tags = entry.tags + [("AL",1)]
            else:
                entry.tags = entry.tags + [("AL",0)]
            outBAM.write(entry)
        outBAM.close()
        inBAM.close()
        return 
    
    if(args.paired_end==True or args.paired_end=="True"):
        if(len(args.fq) != 2):
            raise NameError("I DON'T KNOW WHAT WE'RE YELLING ABOUT. REQUIRED: TWO FASTQ FILES")
        read1 = SeqIO.parse(args.fq[0], "fastq")
        read2 = SeqIO.parse(args.fq[1], "fastq")
        #inBAM = removeSecondary(args.bam_file) #Artefactual code
        postFilterBAM = Samfile(args.bam_file,"rb")
        output=args.output
        if(output=="default"):
            output='.'.join(args.bam_file.split('.')[0:-1])+'.tagged.bam'
        outBAM = Samfile(output,"wb",template=postFilterBAM)
        for entry in postFilterBAM:
            if(entry.is_secondary):
                continue
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
    print("Attempting to remove secondary")
    from subprocess import call
    if(outBAM=="default"):
        outBAM = '.'.join(inBAM.split('.')[0:-1])+'.2ndrm.bam'
    command = "samtools view -hb -F 0x0100 {} -o {}".format(inBAM,outBAM)
    call(command,shell=True)
    print(command + "is shell call")
    return outBAM

if(__name__=="__main__"):
    main()
