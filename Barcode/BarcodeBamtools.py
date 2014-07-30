#!/mounts/anaconda/bin/python

from pysam import Samfile
from Bio import SeqIO
from BarcodeUtils import hamming

def Bam2Sam(inbam,outsam):
    from subprocess import call
    output = open(outsam,'w',0)
    command_str='samtools view -h {}'.format(inbam)
    print(command_str)
    call(command_str,stdout=output,shell=True)
    return(command_str,outsam)

def fuzzyJoining(input,trueFamilyList,output="default"):
    if(output=="default"):
        output = '.'.join(input.split('.')[0:-1]) + ".fuzzy.bam"
    inputBAM=Samfile(input,"rb")
    outputBAM=Samfile(output,"wb",template=inputBAM)
    tagsFile = open(trueFamilyList,"r")
    superFamilyTags = [line.strip().split() for line in tagsFile]
    for entry in inputBAM:
        BSloc = [i for i,j in enumerate(entry.tags) if j[0]=="BS"][0]
        if entry.tags[BSloc][1] in superFamilyTags:
            outputBAM.write(entry)
        else:
            print()
            hamDistances = [hamming(entry.tags[BSloc][1],familyTag) for familyTag in superFamilyTags]
            minDist = min(hamDistances)
            if(minDist<3):
                barcodePos = hamDistances.index(minDist)
                entry.tags[BSloc]=[("BS",superFamilyTags[barcodePos])]
                entry.tags = entry.tags + [("BD",minDist)]
            outputBAM.write(entry)
    return output

def GenerateBarcodeIndexBAM(doubleTaggedBAM,index_file="default"):
    from subprocess import call
    if(index_file=="default"):
        index_file = '.'.join(doubleTaggedBAM.split('.')[0:-2]) + ".DoubleIdx"
    call("samtools view {} | grep -v 'AL:i:0' | awk '{{print $(NF-2)}}' | sed 's/BS:Z://g' | sort | uniq -c | awk 'BEGIN {{OFS=\"\t\"}};{{print $1,$2}}' > {}".format(doubleTaggedBAM,index_file),shell=True)
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
        try:
            famSize = BarDict[record.tags[BSloc][1]]
        except KeyError:
            famSize = 0
        try:
            FMloc = [i for i,j in enumerate(record.tags) if j[0]=="FM"][0]
            record.tags[FMloc] = [("FM",famSize)]
        except IndexError:
            record.tags = record.tags + [("FM",famSize)]
        if(int(famSize) >= 4):
            writeList.write('{}\n'.format(record.tags[BSloc][1]))
            record.tags = record.tags + [("BD",0)]
        outfile.write(record)
    writeList.close()
    outfile.close()
    from subprocess import call
    call("cat {} | sort | uniq > omgz.tmp;mv omgz.tmp {}".format(trueFamilyList, trueFamilyList),shell=True) #TODO: make random name for the temp file!
    return output,trueFamilyList

def mergeBarcodes(reads1,reads2,outfile="default"):
    reader1=Samfile(reads1,"rb")
    reader2=Samfile(reads2,"rb")
    if(outfile=="default"):
        outfile = '.'.join(reads1.split('.')[0:-2])+'merged.bam'
    outSAM = Samfile(outfile,"wb",template=reader1)
    for entry1 in reader1:
        entry2=reader2.next()
        assert entry1.qname==entry2.qname
        BSloc1 = [i for i,j in enumerate(entry1.tags) if j[0]=="BS"][0]
        BSloc2 = [i for i,j in enumerate(entry2.tags) if j[0]=="BS"][0]
        Barcode1, Barcode2 = entry1.tags[BSloc1][1], entry2.tags[BSloc2][1]
        #print("Barcode 1, Barcode 2 are {}, {}, respectively.".format(Barcode1,Barcode2))
        concatBarcode=Barcode1+Barcode2
        #print("New barcode will be {}".format(concatBarcode))
        entry1.tags = entry1.tags[0:BSloc1] + [("BS",concatBarcode)] + entry1.tags[BSloc1+1:]
        entry2.tags = entry2.tags[0:BSloc2] + [("BS",concatBarcode)] + entry2.tags[BSloc2+1:]
        outSAM.write(entry1)
        outSAM.write(entry2)
    reader1.close()
    reader2.close()
    outSAM.close()
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

#Filters out both reads in a pair based on at least one criterion. Both reads must pass to be written.
#Required: SAM file be sorted by name, supplementary and secondary alignments removed, unmapped reads retained.
def pairedFilterBam(inputBAM, passBAM="default", failBAM="default", criteria="default"):
    if(criteria=="default"):
        raise NameError("Filter Failed: Criterion Not Set. Currently supported: adapter pass/fail (\"adapter\"), mapped")
    if(passBAM=="default"):
        passBAM = '.'.join(inputBAM.split('.')[0:-1]) + ".{}.pass.bam".format(criteria)
    if(failBAM=="default"):
        failBAM = '.'.join(inputBAM.split('.')[0:-1]) + ".{}.fail.bam".format(criteria)
    inBAM=Samfile(inputBAM,"rb")
    passFilter = Samfile(passBAM,"wb",template=inBAM)
    failFilter = Samfile(failBAM,"wb",template=inBAM)

    if("adapter" in criteria.lower()):
        while True:
            try:
                read1=inBAM.next()
                read2=inBAM.next() #TODO: add sanity check here to make sure that the reads are each other's mates
                assert read1.qname==read2.qname
                ALloc1 = [i for i,j in enumerate(read1.tags) if j[0]=="AL"][0]
                ALValue1 = int(read2.tags[ALloc1][1])
                ALloc2 = [i for i,j in enumerate(read2.tags) if j[0]=="AL"][0]
                ALValue2 = int(read2.tags[ALloc2][1])
                if(ALValue1==0 or ALValue2==0):
                    failFilter.write(read1)
                    failFilter.write(read2)
                else:
                    passFilter.write(read1)
                    passFilter.write(read2)
            except IndexError:
                break
        passFilter.close()
        failFilter.close()
        inBAM.close()
        return passBAM, failBAM
    
    if("editdistance" in criteria.lower()):
        while True:
            try:
                read1=inBAM.next()
                read2=inBAM.next() #TODO: add sanity check here to make sure that the reads are each other's mates
                assert read1.qname==read2.qname
                NMloc1 = [i for i,j in enumerate(read1.tags) if j[0]=="NM"][0]
                NMloc2 = [i for i,j in enumerate(read2.tags) if j[0]=="NM"][0]
                if(read1.tags[NMloc1][1] == 0 or read2.tags[NMloc2] == 0):    
                    failFilter.write(read1)
                    failFilter.write(read2)
                else:
                    passFilter.write(read1)
                    passFilter.write(read2)
            except IndexError:
                break
        passFilter.close()
        failFilter.close()
        inBAM.close()
        return passBAM, failBAM
    
    if("family" in criteria.lower()):
        while True:
            try:
                read1=inBAM.next()
                read2=inBAM.next() #TODO: add sanity check here to make sure that the reads are each other's mates
                assert read1.qname==read2.qname
                FMloc1 = [i for i,j in enumerate(read1.tags) if j[0]=="FM"][0]
                FMValue1 = int(read1.tags[FMloc1][1])
                if(FMValue1 < 3):
                    failFilter.write(read1)
                    failFilter.write(read2)
                else:
                    passFilter.write(read1)
                    passFilter.write(read2)
            except IndexError:
                break
        passFilter.close()
        failFilter.close()
        inBAM.close()
        return passBAM, failBAM
    
    if("fuzzy" in criteria.lower()):
        while True:
            try:
                read1=inBAM.next()
                read2=inBAM.next() #TODO: add sanity check here to make sure that the reads are each other's mates
                assert read1.qname==read2.qname
                try:
                    BDloc1 = [i for i,j in enumerate(read1.tags) if j[0]=="BD"][0]
                    passFilter.write(read1)
                    passFilter.write(read2)
                except IndexError:
                    failFilter.write(read1)
                    failFilter.write(read2)
            except IndexError:
                break
        passFilter.close()
        failFilter.close()
        inBAM.close()
        return passBAM, failBAM

    
    if("ismapped" in criteria.lower()):
        while True:
            try:
                read1=inBAM.next()
                read2=inBAM.next() #TODO: add sanity check here to make sure that the reads are each other's mates
                assert read1.qname==read2.qname
                if(read1.is_unmapped or read2.is_unmapped):
                    failFilter.write(read1)
                    failFilter.write(read2)
                else:
                    passFilter.write(read1)
                    passFilter.write(read2)
            except IndexError:
                break
        passFilter.close()
        failFilter.close()
        inBAM.close()
        return passBAM, failBAM
 
    if("qc" in criteria.lower()):
        while True:
            try:
                read1=inBAM.next()
                read2=inBAM.next() #TODO: add sanity check here to make sure that the reads are each other's mates
                assert read1.qname==read2.qname
                if(read1.is_qcfail or read2.is_qcfail):
                    failFilter.write(read1)
                    failFilter.write(read2)
                else:
                    passFilter.write(read1)
                    passFilter.write(read2)
            except IndexError:
                break
        passFilter.close()
        failFilter.close()
        inBAM.close()
        return passBAM, failBAM
   
    raise NameError("No valid sorting option selected!")    
    return

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
