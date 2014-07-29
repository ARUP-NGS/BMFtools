#!/mounts/anaconda/bin/python

import pysam
from pysam import Samfile
from Bio import SeqIO
import argparse

def Bam2Sam(inbam,outsam):
    from subprocess import call
    output = open(outsam,'w',0)
    command_str='samtools view -h {}'.format(inbam)
    print(command_str)
    call(command_str,stdout=output,shell=True)
    return(command_str,outsam)

def GenerateBarcodeIndexBAM(doubleTaggedBAM,index_file="default"):
    from subprocess import call
    if(index_file=="default"):
        index_file = '.'.join(doubleTaggedBAM.split('.')[0:-2]) + ".DoubleIdx"
    call("samtools view {} | awk '{print $(NF-2)}' | sed 's/BS:Z://g' | sort | uniq -c | awk 'BEGIN {{OFS=\"\t\"}};{{print $1,$2}}' > {}".format(doubleTaggedBAM,index_file),shell=True)
    return index_file


def getFamilySizeBAM(inputBAM,indexFile,output="default",trueFamilyList="default"):
    if(output=="default"):
        output = inputBAM.split('.')[0] + ".doubleFam.bam"
    if(trueFamilyList=="default"):
        trueFamilyList = inputBAM.split('.')[0] + ".doubleFam.lst"
    index = open(indexFile,"r")
    dictEntries = [line.split() for line in index]
    BarDict = {}
    for entry in dictEntries:
        BarDict[entry[1]]=entry[0]
    index.close()
    Reads=Samfile(inputBAM,"rb")
    outfile = Samfile(output,"wb",template=Reads)
    writeList = open(trueFamilyList,"w",0)
    for record in Reads:
        BSloc = [i for i,j in enumerate(record.tags) if j[0]=="BS"][0]
        famSize = BarDict[record.tags[BSloc][1]]
        try:
            FMloc = [i for i,j in enumerate(record.tags) if j[0]=="FM"][0]
            record.tags[FMloc] = [("FM",famSize)]
        except IndexError:
            record.tags = record.tags + [("FM",famSize)]
        outfile.write(record)
        if(int(famSize) > 2):
            writeList.write('{}\n'.format(record.tags[BSloc][1]))
    writeList.close()
    outfile.close()
    from subprocess import call
    call("cat {} | sort | uniq > omgz.tmp;mv omgz.tmp {}".format(trueFamilyList, trueFamilyList),shell=True) #TODO: make random name for the temp file!
    return output,trueFamilyList

def mergeBarcodes(reads1,reads2,outfile="default",doubleBarcodeList="default"):
    reader1=Samfile(reads1,"rb")
    reader2=Samfile(reads2,"rb")
    if(outfile=="default"):
        outfile = '.'.join(reads1.split('.')[0:-2])+'merged.bam'
    if(doubleBarcodeList=="default"):
        doubleBarcodeList = '.'.join(reads1.split('.')[0:-2])+'merge.lst'
    outSAM = Samfile(outfile,"wb",template=reader1)
    doubleMergedList = open(doubleBarcodeList,"w",0)
    for entry1 in reader1:
        entry2=reader2.next()
        BSloc1 = [i for i,j in enumerate(entry1.tags) if j[0]=="BS"][0]
        BSloc2 = [i for i,j in enumerate(entry2.tags) if j[0]=="BS"][0]
        Barcode1, Barcode2 = entry1.tags[BSloc1][1], entry2.tags[BSloc2][1]
        #print("Barcode 1, Barcode 2 are {}, {}, respectively.".format(Barcode1,Barcode2))
        concatBarcode=Barcode1+Barcode2
        entry1.tags[BSloc1]=["BS",concatBarcode]
        entry2.tags[BSloc2]=["BS",concatBarcode]
        outSAM.write(entry1)
        outSAM.write(entry2)
        doubleMergedList.write(concatBarcode+"\n")
    reader1.close()
    reader2.close()
    outSAM.close()
    doubleMergedList.close()
    return outfile

def pairedBarcodeTagging(fq1,fq2,bam,outputBAM="default"):
    if(outputBAM=="default"):
        outputBAM = '.'.join(bam.split('.')[0:-1]) + "tagged.bam"
    read1 = SeqIO.parse(fq1, "fastq")
    read2 = SeqIO.parse(fq2, "fastq")
    #inBAM = removeSecondary(args.bam_file) #Artefactual code
    postFilterBAM = Samfile(bam,"rb")
    outBAM = Samfile(outputBAM,"wb",template=postFilterBAM)
    for entry in postFilterBAM:
        if(entry.is_secondary or entry.flag > 2048):
            continue
        if(entry.is_read1):
            #print("Read is read 1, with name {}. Now getting variables from read 1's files for the tags.".format(entry.qname))
            tempRead = read1.next()
            #print("Read description (for debugging purposes) is {}".format(tempRead.description))
        elif(entry.is_read2):
            #print("Read is read 2, with name {}. Now getting variables from read 2's files for the tags.".format(entry.qname))
            tempRead = read2.next()
            #print("Read description (for debugging purposes) is {}".format(tempRead.description))
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
    return outputBAM

def removeSecondary(inBAM,outBAM="default"):
    if(outBAM=="default"):
        outBAM = '.'.join(inBAM.split('.')[0:-1])+'.2ndrm.bam'
    input=Samfile(inBAM,"rb")
    output=Samfile(outBAM,"wb",template=input)
    print("Attempting to remove secondary")
    for entry in input:
        if(entry.is_secondary or entry.flag > 2048):
            continue
        else:
            output.write(entry)
    return outBAM

def Sam2Bam(insam,outbam):
    from subprocess import call
    output = open(outbam,'w',0)
    command_str='samtools view -Sbh {}'.format(insam,shell=True)
    print(command_str)
    call(command_str,stdout=output,shell=True)
    return(command_str,outbam)

def singleBarcodeTagging(fastq,bam,outputBAM="default"):
    if(outputBAM=="default"):
        outputBAM = '.'.join(bam.split('.')[0:-1]) + "tagged.bam"
    print("Tagged BAM file is: {}.".format(outputBAM))
    reads = SeqIO.parse(fastq, "fastq")
    #inBAM = removeSecondary(args.bam_file) #Artefactual code
    postFilterBAM = Samfile(bam,"rb")
    outBAM = Samfile(outputBAM,"wb",template=postFilterBAM)
    for entry in postFilterBAM:
        if(entry.is_secondary or entry.flag>2048):
            continue
        else:
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
    postFilterBAM.close()
    return outputBAM

def splitBAMByReads(BAM,read1BAM="default",read2BAM="default"):
    presplitBAM = Samfile(BAM,"rb")
    if(read1BAM=="default"):
        read1BAM = '.'.join(BAM.split('.')[0:-1])+'.R1.bam'
    if(read2BAM=="default"):
        read2BAM = '.'.join(BAM.split('.')[0:-1])+'.R2.bam'
    out1 = Samfile(read1BAM,"wb",template=presplitBAM)
    out2 = Samfile(read2BAM,"wb",template=presplitBAM)
    for entry in presplitBAM:
        if(entry.is_read1):
            out1.write(entry)
        if(entry.is_read2):
            out2.write(entry)
    out1.close()
    out2.close()
    return read1BAM, read2BAM

def pairedFilterFamilyBam(inputBAM,failedBAM):
    
    return