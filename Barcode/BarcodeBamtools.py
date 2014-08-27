#!/mounts/anaconda/bin/python

import pysam
from Bio import SeqIO


def Bam2Sam(inbam, outsam):
    from subprocess import call
    output = open(outsam, 'w', 0)
    command_str = 'samtools view -h {}'.format(inbam)
    print(command_str)
    call(command_str, stdout=output, shell=True)
    return(command_str, outsam)

def BarcodeSort(inbam, outbam="default",paired=True):
    if(outbam == "default"):
        outbam = '.'.join(inbam.split('.')[0:-1]) + "barcodeSorted.bam"
    outsam = '.'.join(outbam.split('.')[0:-1]) + "barcodeSorted.sam"
    from subprocess import call
    call("samtools view -H {} > {}".format(inbam, outsam), shell=True)
    print("Now converting bam to sam for sorting by barcode.")
    if(paired==False):
        call("samtools view {} | awk 'BEGIN {{FS=\"\t\";OFS=\"\t\"}};{{print $(NF-2),$0}}' - | sort | cut -f2- -d' ' >> {}".format(inbam, outsam), shell=True)
    else:
        call("samtools view {} | awk 'BEGIN {{FS=\"\t\";OFS=\"\t\"}};{{print $(NF-1),$0}}' - | sort | cut -f2- -d' ' >> {}".format(inbam, outsam), shell=True)
    print("Now converting sam back to bam for further operations.")
    call("samtools view -Sbh {} > {}".format(outsam, outbam), shell=True)
    call("rm {}".format(outsam), shell=True)
    return outbam

def compareSamRecords(RecordList):
    seqs = [record.seq for record in RecordList]
    max = 0
    Success = False
    for seq in seqs:
        numEq = sum(seq == seqItem for seqItem in seqs)
        if(numEq > max):
            max = numEq
            finalSeq = seq
    frac = numEq * 1.0 / len(RecordList)
    if(frac > 0.9):
        Success = True
    consolidatedRecord = RecordList[0]
    consolidatedRecord.seq = finalSeq
    return consolidatedRecord, Success

def Consolidate(inbam, outbam="default"):
    if(outbam == "default"):
        outbam = '.'.join(inbam.split('.')[0:-1]) + "consolidated.bam"
    inputHandle = pysam.Samfile(inbam, 'rb')
    outputHandle = pysam.Samfile(outbam, 'wb', template=inputHandle)
    workingBarcode1 = ""
    workingBarcode2 = ""
    workingSet1 = []
    workingSet2 = []
    for record in inputHandle:
        if(record.is_read1):
            BSloc1 = [i for i, j in enumerate(record.tags) if j[0] == "BS"][0]
            barcodeRecord1 = record.tags[BSloc1][1]
            if(workingBarcode1 == ""):
                workingBarcode1 = barcodeRecord1
                workingSet1 = []
                workingSet1.append(record)
            elif(workingBarcode1 == barcodeRecord1):
                workingSet1.append(record)
            else:
                mergedRecord1, success = compareSamRecords(workingSet1)
                if(success == True):
                    outputHandle.write(mergedRecord1)
                WorkingSet1 = []
                WorkingSet1.append(record)
                workingBarcode1 = barcodeRecord1
        if(record.is_read2):
            BSloc2 = [i for i, j in enumerate(record.tags) if j[0] == "BS"][0]
            barcodeRecord2 = record.tags[BSloc2][1]
            if(workingBarcode2 == ""):
                workingBarcode2 = barcodeRecord2
                workingSet2 = []
                workingSet2.append(record)
            elif(workingBarcode2 == barcodeRecord2):
                workingSet2.append(record)
            else:
                mergedRecord2, success = compareSamRecords(workingSet2)
                if(success == True):
                    outputHandle.write(mergedRecord2)
                WorkingSet2 = []
                WorkingSet2.append(record)
                workingBarcode2 = barcodeRecord2
    inputHandle.close()
    outputHandle.close()
    return outbam 

def CorrSort(inbam, outprefix="default"):
    from subprocess import call
    print("inbam variable is {}".format(inbam))
    if(outprefix == "default"):
        outprefix = '.'.join(inbam.split('.')[0:-1]) + ".CorrSort"
    command_str = "samtools sort {} {}".format(inbam, outprefix)
    print(command_str)
    call(command_str, shell=True)
    return(outprefix + ".bam")

def criteriaTest(read1, read2, filter="default"):
    list = "adapter barcode complexity editdistance family ismapped qc".split(' ')
    
    if(filter == "default"):
        print("List of valid filters: {}".format(', '.join(list)))
        raise ValueError("Filter must be set! Requires an exact match (case insensitive).")
    
    if(filter not in list):
        raise ValueError("Your filter is not a supported criterion. The list is: {}".format(list))
    
    if(filter == "adapter"):
        ALloc1 = [i for i, j in enumerate(read1.tags) if j[0] == "AL"][0]
        try:
            ALValue1 = int(read1.tags[ALloc1][1])
        except IndexError:
            ALValue1 = 0
        ALloc2 = [i for i, j in enumerate(read2.tags) if j[0] == "AL"][0]
        try:
            ALValue2 = int(read2.tags[ALloc2][1])
        except IndexError:
            ALValue2 = 0
        if(ALValue1 == 0 or ALValue2 == 0 or ALValue1 == "0" or ALValue2 == "0"):
            return False
    
    if(filter == "barcode"):
        BSloc1 = [i for i, j in enumerate(read1.tags) if j[0] == "BS"][0] 
        BSloc2 = [i for i, j in enumerate(read2.tags) if j[0] == "BS"][0]
        if(read1.tags[BSloc1][1] != read2.tags[BSloc2][1]):
            return False
        
    if(filter == "complexity"):
        if("AAAAAAAAAAA" in str(read1.tags) + str(read2.tags) or "TTTTTTTTTTT" in str(read1.tags) + str(read2.tags) or "GGGGGGGGGGG" in str(read1.tags) + str(read2.tags) or "CCCCCCCCCCC" in str(read1.tags) + str(read2.tags)):
            return False
        
    if(filter == "editdistance"):
        NMloc1 = [i for i, j in enumerate(read1.tags) if j[0] == "NM"][0]
        NMloc2 = [i for i, j in enumerate(read2.tags) if j[0] == "NM"][0]
        if(read1.tags[NMloc1][1] == 0 and read2.tags[NMloc2] == 0):
            return False
        
    if(filter == "family"):
        FMloc1 = [i for i, j in enumerate(read1.tags) if j[0] == "FM"][0]
        FMloc2 = [i for i, j in enumerate(read2.tags) if j[0] == "FM"][0]
        FMValue1 = int(read1.tags[FMloc1][1])
        FMValue2 = int(read2.tags[FMloc2][1])
        if(FMValue1 < 4 or FMValue2 < 4):
            return False
    
    if(filter == "ismapped"):
        if(read1.is_unmapped or read2.is_unmapped):
            return False
     
    if(filter == "qc"):
        if(read1.is_qcfail or read2.is_qcfail):
            return False
    
    return True
        
# I'm not a huge fan of this. We could only rescue a small fraction. Maybe it becomes valuable if we have deeper coverage, but for now, we can't really get much out of this.
'''
def fuzzyJoining(inputBAM,referenceFasta,output="default",):
    import subprocess
    return
'''

def GenerateBarcodeIndexBAM(doubleTaggedBAM, index_file="default",paired=True):
    from subprocess import call
    if(index_file == "default"):
        index_file = '.'.join(doubleTaggedBAM.split('.')[0:1]) + ".DoubleIdx"
    if(paired==True):
        call("samtools view {} | grep -v 'AL:i:0' | awk '{{print $(NF-1)}}' | sed 's/BS:Z://g' | sort | uniq -c | awk 'BEGIN {{OFS=\"\t\"}};{{print $1,$2}}' > {}".format(doubleTaggedBAM, index_file), shell=True)
    if(paired==False):
        call("samtools view {} | grep -v 'AL:i:0' | awk '{{print $(NF-2)}}' | sed 's/BS:Z://g' | sort | uniq -c | awk 'BEGIN {{OFS=\"\t\"}};{{print $1,$2}}' > {}".format(doubleTaggedBAM, index_file), shell=True)
    return index_file

def GenerateBarcodeIndexReference(doubleIdx, output="default"):
    if(output == "default"):
        output = doubleIdx.split('.')[0] + '.ref.fasta'
    from subprocess import call
    call("paste {0} {0} | sed 's:^:>:g' | tr '\t' '\n' > {1}".format(doubleIdx, output), shell=True)
    return output

def GenerateFamilyHistochart(BamBarcodeIndex, output="default"):
    from subprocess import call
    if(output == "default"):
        output = '.'.join(BamBarcodeIndex.split('.')[:-1]) + 'hist.txt'
    commandStr = "cat {} | awk '{{print $1}}' | sort | uniq -c | awk 'BEGIN {{OFS=\"\t\"}};{{print $1,$2}}' >> {}".format(BamBarcodeIndex, output)
    print("Command str: {}".format(commandStr))
    call(commandStr, shell=True)
    print("Histochart for family sizes now available at {}".format(output))
    return output

def getFamilySizeBAM(inputBAM, indexFile, output="default", trueFamilyList="default"):
    if(output == "default"):
        output = inputBAM.split('.')[0] + ".doubleFam.bam"
    if(trueFamilyList == "default"):
        trueFamilyList = inputBAM.split('.')[0] + ".doubleFam.lst"
    index = open(indexFile, "r")
    dictEntries = [line.split() for line in index]
    BarDict = {}
    for entry in dictEntries:
        BarDict[entry[1]] = entry[0]
    index.close()
    Reads = pysam.Samfile(inputBAM, "rb")
    outfile = pysam.Samfile(output, "wb", template=Reads)
    writeList = open(trueFamilyList, "w", 0)
    for record in Reads:
        BSloc = [i for i, j in enumerate(record.tags) if j[0] == "BS"][0]
        try:
            famSize = BarDict[record.tags[BSloc][1]]
        except KeyError:
            famSize = 0
        try:
            FMloc = [i for i, j in enumerate(record.tags) if j[0] == "FM"][0]
            record.tags = record.tags[0:FMloc] + [("FM", famSize)] + record.tags[FMloc + 1:]
        except IndexError:
            record.tags = record.tags + [("FM", famSize)]
        if(int(famSize) >= 2):
            writeList.write('{}\n'.format(record.tags[BSloc][1]))
            try:
                BDloc = [i for i, j in enumerate(record.tags) if j[0] == "BD"][0]
                if(BDloc == len(record.tags) - 1):
                    record.tags = record.tags[0:BDloc] + [("BD", 0)] 
                else:
                    record.tags = record.tags[0:BDloc] + [("BD", 0)] + record.tags[BDloc + 1:]
            except IndexError:
                record.tags = record.tags + [("BD", 0)]
        outfile.write(record)
    writeList.close()
    outfile.close()
    from subprocess import call
    import uuid
    tempname = str(uuid.uuid4().get_hex().upper()[0:12]) + ".OMGZZZZ.tmp"
    call("cat {0} | sort | uniq > {1};mv {1} {0}".format(trueFamilyList, tempname), shell=True)
    return output, trueFamilyList

def mergeBarcodes(reads1, reads2, outfile="default"):
    reader1 = pysam.Samfile(reads1, "rb")
    reader2 = pysam.Samfile(reads2, "rb")
    if(outfile == "default"):
        outfile = '.'.join(reads1.split('.')[0:-2]) + '.merged.bam'
    outSAM = pysam.Samfile(outfile, "wb", template=reader1)
    for entry1 in reader1:
        entry2 = reader2.next()
        assert entry1.qname == entry2.qname
        BSloc1 = [i for i, j in enumerate(entry1.tags) if j[0] == "BS"][0]
        BSloc2 = [i for i, j in enumerate(entry2.tags) if j[0] == "BS"][0]
        Barcode1, Barcode2 = entry1.tags[BSloc1][1], entry2.tags[BSloc2][1]
        # print("Barcode 1, Barcode 2 are {}, {}, respectively.".format(Barcode1,Barcode2))
        concatBarcode = Barcode1 + Barcode2
        # print("New barcode will be {}".format(concatBarcode))
        entry1.tags = entry1.tags[0:BSloc1] + [("BS", concatBarcode)] + entry1.tags[BSloc1 + 1:]
        entry2.tags = entry2.tags[0:BSloc2] + [("BS", concatBarcode)] + entry2.tags[BSloc2 + 1:]
        outSAM.write(entry1)
        outSAM.write(entry2)
    reader1.close()
    reader2.close()
    outSAM.close()
    return outfile

def pairedBarcodeTagging(fq1, fq2, bam, outputBAM="default",secondSuppBAM="default"):
    if(outputBAM == "default"):
        outputBAM = '.'.join(bam.split('.')[0:-1]) + "tagged.bam"
    if(secondSuppBAM=="default"):
        secondSuppBAM = bam.split('.')[0]+'.2ndSupp.bam'
    read1 = SeqIO.parse(fq1, "fastq")
    read2 = SeqIO.parse(fq2, "fastq")
    # inBAM = removeSecondary(args.bam_file) #Artefactual code
    postFilterBAM = pysam.Samfile(bam, "rb")
    outBAM = pysam.Samfile(outputBAM, "wb", template=postFilterBAM)
    suppBAM = pysam.Samfile(secondSuppBAM,"wb",template=postFilterBAM)
    for entry in postFilterBAM:
        if(entry.is_secondary or entry.flag > 2048):
            suppBAM.write(entry)
            continue
        if(entry.is_read1):
            # print("Read is read 1, with name {}. Now getting variables from read 1's files for the tags.".format(entry.qname))
            tempRead = read1.next()
            # print("Read description (for debugging purposes) is {}".format(tempRead.description))
        elif(entry.is_read2):
            # print("Read is read 2, with name {}. Now getting variables from read 2's files for the tags.".format(entry.qname))
            tempRead = read2.next()
        descArray = tempRead.description.split("###")
        try:
            entry.tags = entry.tags + [("BS", descArray[2].strip())]
        except IndexError:
            print(descArray)
            raise ValueError("Something's wrong!!!")
        if(descArray[1].strip() == "AdapterPass"):
            entry.tags = entry.tags + [("AL", 1)]
        else:
            entry.tags = entry.tags + [("AL", 0)]
        outBAM.write(entry)
    suppBAM.close()
    outBAM.close()
    postFilterBAM.close()
    return outputBAM

# Filters out both reads in a pair based on a list of comma-separated criteria. Both reads must pass to be written.
# Required: SAM file be sorted by name, supplementary and secondary alignments removed, unmapped reads retained.
def pairedFilterBam(inputBAM, passBAM="default", failBAM="default", criteria="default"):
    if(criteria == "default"):
        raise NameError("Filter Failed: Criterion Not Set. Currently supported: adapter pass/fail (\"adapter\"),\"ismapped\",\"editdistance\",\"family\",\"fuzzy\",\"qc\"")
    if(passBAM == "default"):
        passBAM = '.'.join(inputBAM.split('.')[0:-1]) + ".{}P.bam".format(criteria)
    if(failBAM == "default"):
        failBAM = '.'.join(inputBAM.split('.')[0:-1]) + ".{}F.bam".format(criteria)
    inBAM = pysam.Samfile(inputBAM, "rb")
    passFilter = pysam.Samfile(passBAM, "wb", template=inBAM)
    failFilter = pysam.Samfile(failBAM, "wb", template=inBAM)
    criteriaList = criteria.lower().split(',')
    for i, entry in enumerate(criteriaList):
        print("Criteria #{} is \"{}\"".format(i, entry))
    for read in inBAM:
        failed = False
        if(read.is_read1):
            read1 = read
            continue
        if(read.is_read2):
            read2 = read
            assert read1.qname == read2.qname  # Sanity check here to make sure that the reads are each other's mates
            for criterion in criteriaList:
                result = criteriaTest(read1, read2, filter=criterion)
                if(result == False):
                    failed = True
                    failFilter.write(read1)
                    failFilter.write(read2)
                    break
                else:
                    continue
            if(failed == True):
                continue
            passFilter.write(read1)
            passFilter.write(read2)
    passFilter.close()
    failFilter.close()
    inBAM.close()
    return passBAM, failBAM

def removeSecondary(inBAM, outBAM="default"):
    if(outBAM == "default"):
        outBAM = '.'.join(inBAM.split('.')[0:-1]) + '.2ndrm.bam'
    input = pysam.Samfile(inBAM, "rb")
    output = pysam.Samfile(outBAM, "wb", template=input)
    print("Attempting to remove secondary")
    for entry in input:
        if(entry.is_secondary or entry.flag > 2048):
            continue
        else:
            output.write(entry)
    return outBAM

def Sam2Bam(insam, outbam):
    from subprocess import call
    output = open(outbam, 'w', 0)
    command_str = 'samtools view -Sbh {}'.format(insam, shell=True)
    print(command_str)
    call(command_str, stdout=output, shell=True)
    return(command_str, outbam)

# Taken from Scrutils, by Jacob Durtschi
def SamtoolsBam2fq(bamPath, outFastqPath):
    import subprocess
    # Build commands that will be piped
    samtoolsBam2fqCommand = [ 
                'samtools', 'bam2fq',
                bamPath
                ]   

    gzipCommand = [ 
                'gzip'
                ]   

    # Open output fastq.gz file to be piped into
    outFastq = open(outFastqPath, 'w')

    # Call piped commands
    process1 = subprocess.Popen(samtoolsBam2fqCommand, stdout=subprocess.PIPE, shell=False)
    process2 = subprocess.Popen(gzipCommand, stdin=process1.stdout, stdout=outFastq, shell=False)
    # process1.stdout.close()
    # process1.wait()
    process2.wait()
    outFastq.flush()
    outFastq.close()
    # if processErr:
    #    print('Error: in bwa mem function', end='\n', file=sys.stderr)
    #    print('From bwa and samtools:', end='\n', file=sys.stderr)
    #    print(processErr, end='\n', file=sys.stderr)
    #    sys.exit(1)

    return outFastqPath

def singleBarcodeTagging(fastq, bam, outputBAM="default"):
    if(outputBAM == "default"):
        outputBAM = '.'.join(bam.split('.')[0:-1]) + "tagged.bam"
    print("Tagged BAM file is: {}.".format(outputBAM))
    reads = SeqIO.parse(fastq, "fastq")
    # inBAM = removeSecondary(args.bam_file) #Artefactual code
    postFilterBAM = pysam.Samfile(bam, "rb")
    outBAM = pysam.Samfile(outputBAM, "wb", template=postFilterBAM)
    for entry in postFilterBAM:
        if(entry.is_secondary or entry.flag > 2048):
            continue
        else:
            try:
                tempRead = reads.next()
            except StopIteration:
                break
        descArray = tempRead.description.split("###")
        entry.tags = entry.tags + [("BS", descArray[-2].strip())]
        entry.tags = entry.tags + [("FM", descArray[-1].strip())]
        if("Adapter" not in descArray[-3]):
            print("The value in descArray[-3] is not what it should be! Value: {}".format(descArray[-3]))
            print("descArray is {}".format(descArray))
            raise ValueError("Something has gone wrong! The adapter pass/fail ")
        if(descArray[-3].strip() == "AdapterPass"):
            entry.tags = entry.tags + [("AL", 1)]
        else:
            entry.tags = entry.tags + [("AL", 0)]
        outBAM.write(entry)
    outBAM.close()
    postFilterBAM.close()
    return outputBAM


def SingleConsolidate(inbam, outbam="default"):
    if(outbam == "default"):
        outbam = '.'.join(inbam.split('.')[0:-1]) + "consolidated.bam"
    inputHandle = pysam.Samfile(inbam, 'rb')
    outputHandle = pysam.Samfile(outbam, 'wb', template=inputHandle)
    workingBarcode = ""
    workingSet = []
    for record in inputHandle:
        BSloc = [i for i, j in enumerate(record.tags) if j[0] == "BS"][0]
        barcodeRecord = record.tags[BSloc][1]
        if(workingBarcode == ""):
            workingBarcode = barcodeRecord
            workingSet = []
            workingSet.append(record)
        elif(workingBarcode == barcodeRecord):
            workingSet.append(record)
        else:
            mergedRecord, success = compareSamRecords(workingSet)
            if(success == True):
                outputHandle.write(mergedRecord)
            WorkingSet = []
            WorkingSet.append(record)
            workingBarcode = barcodeRecord
    inputHandle.close()
    outputHandle.close()
    return outbam 

def singleCriteriaTest(read, filter="default"):
    list = "adapter complexity editdistance family ismapped qc".split(' ')
    
    if(filter == "default"):
        print("List of valid filters: {}".format(', '.join(list)))
        raise ValueError("Filter must be set! Requires an exact match (case insensitive).")
    
    if(filter not in list):
        raise ValueError("Your filter is not a supported criterion. The list is: {}".format(list))
    
    if(filter == "adapter"):
        ALloc1 = [i for i, j in enumerate(read.tags) if j[0] == "AL"][0]
        try:
            ALValue1 = int(read.tags[ALloc1][1])
        except IndexError:
            ALValue1 = 0
        if(ALValue1 == 0):
            return False
    
    if(filter == "complexity"):
        if("AAAAAAAAAAA" in str(read.tags) or "TTTTTTTTTTT" in str(read.tags) or "GGGGGGGGGGG" in str(read.tags) or "CCCCCCCCCCC" in str(read.tags)):
            return False
        
    if(filter == "editdistance"):
        NMloc1 = [i for i, j in enumerate(read.tags) if j[0] == "NM"][0]
        if(read.tags[NMloc1][1] == 0):
            return False
        
    if(filter == "family"):
        FMloc1 = [i for i, j in enumerate(read.tags) if j[0] == "FM"][0]
        FMValue1 = int(read.tags[FMloc1][1])
        if(FMValue1 < 3):
            return False
    
    if(filter == "ismapped"):
        if(read.is_unmapped):
            return False
     
    if(filter == "qc"):
        if(read.is_qcfail):
            return False
    
    return True

# Filters out both reads in a pair based on a list of comma-separated criteria. Both reads must pass to be written.
# Required: SAM file be sorted by name, supplementary and secondary alignments removed, unmapped reads retained.
def singleFilterBam(inputBAM, passBAM="default", failBAM="default", criteria="default"):
    if(criteria == "default"):
        raise NameError("Filter Failed: Criterion Not Set. Currently supported: adapter pass/fail (\"adapter\"),\"ismapped\",\"editdistance\",\"family\",\"fuzzy\",\"qc\"")
    if(passBAM == "default"):
        passBAM = '.'.join(inputBAM.split('.')[0:-1]) + ".{}P.bam".format(criteria)
    if(failBAM == "default"):
        failBAM = '.'.join(inputBAM.split('.')[0:-1]) + ".{}F.bam".format(criteria)
    inBAM = pysam.Samfile(inputBAM, "rb")
    passFilter = pysam.Samfile(passBAM, "wb", template=inBAM)
    failFilter = pysam.Samfile(failBAM, "wb", template=inBAM)
    criteriaList = criteria.lower().split(',')
    for i, entry in enumerate(criteriaList):
        print("Criteria #{} is \"{}\"".format(i, entry))
    for read in inBAM:
        failed = False
        for criterion in criteriaList:
            result = singleCriteriaTest(read, filter=criterion)
            if(result == False):
                failFilter.write(read)
                failed = True
                break
            else:
                continue
        if(failed == True):
            continue
        passFilter.write(read)
    passFilter.close()
    failFilter.close()
    inBAM.close()
    return passBAM, failBAM

def splitBAMByReads(BAM, read1BAM="default", read2BAM="default"):
    presplitBAM = pysam.Samfile(BAM, "rb")
    if(read1BAM == "default"):
        read1BAM = '.'.join(BAM.split('.')[0:-1]) + '.R1.bam'
    if(read2BAM == "default"):
        read2BAM = '.'.join(BAM.split('.')[0:-1]) + '.R2.bam'
    out1 = pysam.Samfile(read1BAM, "wb", template=presplitBAM)
    out2 = pysam.Samfile(read2BAM, "wb", template=presplitBAM)
    for entry in presplitBAM:
        if(entry.is_read1):
            out1.write(entry)
        if(entry.is_read2):
            out2.write(entry)
    out1.close()
    out2.close()
    return read1BAM, read2BAM
