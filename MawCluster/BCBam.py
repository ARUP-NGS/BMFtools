import shlex
import subprocess
import os
import shutil
import logging

from Bio import SeqIO
import pysam

from MawCluster import BCFastq
from utilBMF.HTSUtils import ThisIsMadness, printlog as pl
from utilBMF import HTSUtils


def AbraCadabra(inbam,
                outbam="default",
                jar="default",
                memStr="default",
                ref="default",
                threads="4",
                bed="default",
                working="default",
                log="default"):
    """
    Calls abra for indel realignment. It supposedly
    out-performs GATK's IndelRealigner.
    Note: bed file must be first 3 columns only and
    coordinate sorted. You will likely need an additional bed file for this.
    """
    if(jar == "default"):
        HTSUtils.FacePalm("Required: Path to abra jar!")
    else:
        pl("Non-default abra jar used: " + jar)
    if(memStr == "default"):
        memStr = "-Xmx16G"
        pl("Default memory string used: " + memStr)
    else:
        pl("Non-default memory string used: " + memStr)
    if(ref == "default"):
        raise ValueError("Reference fasta must be provided!")
    else:
        pl("Reference file set: {}.".format(ref))
    if(bed == "default"):
        raise ValueError("Bed file required.")
    else:
        pl("Bed file set: {}.".format(bed))
    if(working == "default"):
        working = inbam.split('.')[0] + "working_dir"
        pl("Default working directory set to be: " + working)
    else:
        pl("Non-default working directory: " + working)
    if(log == "default"):
        log = "abra.log"
    if(outbam == "default"):
        outbam = '.'.join(inbam.split('.')[0:-1]) + '.abra.bam'
    from os import path
    pl(("Command to reproduce the call of this function: "
        "AbraCadabra({}, outbam={}, jar={}, ".format(inbam, outbam, jar) +
        "memStr={}, ref={}, threads={}, ".format(memStr, ref, threads) +
        "bed={}, working={}, log={})".format(bed, working, log)))
    if(path.isdir(working)):
        pl("Working directory already exists - deleting!")
        shutil.rmtree(working)
    # Check bed file to make sure it is in appropriate format for abra
    bedLines = [line.strip() for line in open(bed, "r").readlines()]
    if(len(bedLines[0]) > 3):
        newbed = '.'.join(bed.split('.')[0:-1]) + 'abra.bed'
        pl("Bed file provided not in form abra accepts.")
        subprocess.check_call("cut -f1-3 {} | sort -k1,1 -k2,2n > {}".format(
                              bed, newbed),
                              shell=True)
        bed = newbed
    import os.path
    if(os.path.isfile(inbam + ".bai") is False):
        pl("No bam index found for input bam - attempting to create.")
        subprocess.check_call(['samtools', 'index', inbam])
        if(os.path.isfile(inbam + ".bai") is False):
            inbam = HTSUtils.CoorSortAndIndexBam(inbam, outbam, uuid=True)
    command = ("java {} -jar {} --in {}".format(memStr, jar, inbam) +
               " --out {} --ref {} --targets".format(outbam, ref) +
               " {} --threads {} ".format(bed, threads) +
               "--working {}".format(working))
    pl("Command: {}.".format(command))
    subprocess.check_call(shlex.split(command), shell=False)
    return outbam


def Bam2Sam(inbam, outsam):
    pl("Bam2Sam. Input: {}. Output: {}.".format(inbam, outsam))
    output = open(outsam, 'w', 0)
    command_str = 'samtools view -h {}'.format(inbam)
    pl(command_str)
    subprocess.check_call(shlex.split(command_str), stdout=output, shell=False)
    return(command_str, outsam)


def BarcodeSort(inbam, outbam="default", paired=True):
    if(outbam == "default"):
        outbam = '.'.join(inbam.split('.')[0:-1]) + "barcodeSorted.bam"
    pl("BarcodeSort. Input: {}. Output: {}.".format(inbam, outbam))
    outsam = '.'.join(outbam.split('.')[0:-1]) + ".sam"
    headerCommand = "samtools view -H {}".format(inbam)
    pl(headerCommand)
    subprocess.check_call(shlex.split(headerCommand),
                          shell=False, stdout=outsam)
    pl("Now converting bam to sam for sorting by barcode.")
    if(paired is False):
        cmd = ("samtools view {} | ".format(inbam) +
               "awk 'BEGIN {{FS=\"\t\";OFS=\"\t\"}};{{print "
               "$(NF-2),$0}}' - | sort | cut -f2- >> {}".format(outsam))
        pl(
            "Command string for this sorting process is: {}".format(cmd))
        HTSUtils.PipedShellCall(cmd)
    else:
        cmd = ("samtools view {} | ".format(inbam) +
               "awk 'BEGIN {{FS=\"\t\";OFS=\"\t\"}};{{print "
               "$(NF-1),$0}}' - | sort | cut -f2- >> {}".format(outsam))
        pl(
            "Command string for this sorting process is: {}".format(cmd))
        HTSUtils.PipedShellCall(cmd)
    pl("Now converting sam back to bam for further operations.")
    subprocess.check_call(shlex.split("samtools view -Sbh {}".format(outsam)),
                          shell=False, stdout=outbam)
    os.remove(outsam)
    return outbam


def compareRecs(RecordList, stringency=0.9):
    seqs = [record.seq for record in RecordList]
    max = 0
    Success = False
    for seq in seqs:
        pl("Seq: {}".format(seq))
        numEq = sum(seq == seqItem for seqItem in seqs)
        if(numEq > max):
            max = numEq
            finalSeq = seq
    frac = numEq * 1.0 / len(RecordList)
    pl("Fraction {}. Stringency: {}. Pass? {}.".format(
        frac, stringency, (frac > stringency)))
    if(frac > stringency):
        Success = True
    consolidatedRecord = RecordList[0]
    consolidatedRecord.seq = finalSeq
    return consolidatedRecord, Success


def Consolidate(inbam, outbam="default", stringency=0.9):
    if(outbam == "default"):
        outbam = '.'.join(inbam.split('.')[0:-1]) + "consolidated.bam"
    inputHandle = pysam.Samfile(inbam, 'rb')
    outputHandle = pysam.Samfile(outbam, 'wb', template=inputHandle)
    workBC1 = ""
    workBC2 = ""
    Set1 = []
    Set2 = []
    for record in inputHandle:
        if(record.is_read1):
            barcodeRecord1 = record.opt("BS")
            if(workBC1 == ""):
                workBC1 = barcodeRecord1
                Set1 = []
                Set1.append(record)
            elif(workBC1 == barcodeRecord1):
                Set1.append(record)
            else:
                mergeRec1, success = compareRecs(Set1, stringency=stringency)
                if(success is True):
                    outputHandle.write(mergeRec1)
                WorkingSet1 = []
                WorkingSet1.append(record)
                workBC1 = barcodeRecord1
        if(record.is_read2):
            barcodeRecord2 = record.opt("BS")
            if(workBC2 == ""):
                workBC2 = barcodeRecord2
                Set2 = []
                Set2.append(record)
            elif(workBC2 == barcodeRecord2):
                Set2.append(record)
            else:
                mergeRec2, success = compareRecs(Set2, stringency=stringency)
                if(success is True):
                    outputHandle.write(mergeRec2)
                WorkingSet2 = []
                WorkingSet2.append(record)
                workBC2 = barcodeRecord2
    inputHandle.close()
    outputHandle.close()
    return outbam


def CoorSort(inbam, outbam="default"):
    import uuid
    pl("CorrSort. Input: {}".format(inbam))
    pl("inbam variable is {}".format(inbam))
    if(outbam == "default"):
        outbam = '.'.join(inbam.split('.')[0:-1]) + ".CoorSort.bam"
    command_str = ("samtools sort {} -o {} -T test".format(inbam, outbam) +
                   str(uuid.uuid4().get_hex().upper()[0:8]))
    pl(command_str)
    subprocess.check_call(shlex.split(command_str), shell=False)
    return(outbam)


def criteriaTest(read1, read2, filterSet="default", minFamSize=3):
    """
    Tool for filtering a pair of BAM files by criteria.
    Note: complexity filter is not needed for shades protocol.
    """
    list = ("adapter barcode complexity editdistance family "
            "ismapped qc notinbed").split()
    Logger = logging.getLogger("Primarylogger")
    try:
        assert isinstance(filterSet, str)  # This should be a string.
    except AssertionError:
        if(isinstance(filterSet, list)):
            pl("filterSet is a list, which was not expected. repr: {}".format(
                repr(filterSet)))
            pl("Joining filterSet into a string.")
            filterSet = ''.join(filterSet)
            if(isinstance(filterSet, str) is False):
                HTSUtils.FacePalm("filterSet isn't a string after joining.")
        else:
            pl("Filter set is not as expected! Repr: {}".format(
                repr(filterSet)))
            raise ValueError("What is this? filterSet is not string or list.")
    if(filterSet == "default"):
        pl("List of valid filters: {}".format(', '.join(list)))
        raise ValueError("Filter must be set!")
    filterSet = filterSet.split(',')
    for filt in filterSet:
        if filt not in list:
            pl("filterSet provided: {}".format(filterSet))
            raise ValueError("Select valid filter(s) - {}".format(list))

    if("adapter" in filterSet):
        ALValue1 = int(read1.opt("FP"))
        ALValue2 = int(read2.opt("FP"))
        if(sum([ALValue1, ALValue2]) != 2):
            return False

    if("barcode" in filterSet):
        if(read1.opt("BS") != read2.opt("BS")):
            return False
            Logger.debug(
                "Barcode sequence didn't match. Are you running shades?")

    if("editdistance" in filterSet):
        NMValue1 = int(read1.opt("NM"))
        NMValue2 = int(read2.opt("NM"))
        if(NMValue1 == 0 and NMValue2 == 0):
            return False

    if("family" in filterSet):
        FMValue1 = int(read1.opt("FM"))
        FMValue2 = int(read2.opt("FM"))
        if(FMValue1 < int(minFamSize) or FMValue2 < int(minFamSize)):
            Logger.debug(("This family didn't survive. Their FMValues:" +
                          "{}, {}".format(FMValue1, FMValue2)))
            return False

    if("ismapped" in filterSet):
        if(read1.is_unmapped and read2.is_unmapped):
            return False

    if("qc" in filterSet):
        if(read1.is_qcfail or read2.is_qcfail):
            return False

    return True


def GenBCIndexBAM(tagBAM, idx="default", paired=True):
    pl("GenBCIndexBAM for {}".format(tagBAM))
    if(idx == "default"):
        idx = '.'.join(tagBAM.split('.')[0:1]) + ".DoubleIdx"
    if(paired is True):
        Str = ("samtools view {} | grep -v 'AL:i:0' | awk ".format(tagBAM) +
               "'{{print $(NF-1)}}' | sed 's/BS:Z://g' | sort | uniq -c |"
               " awk 'BEGIN {{OFS=\"\t\"}};{{print $1,$2}}' > {}".format(idx))
        HTSUtils.PipedShellCall(Str)
    if(paired is False):
        Str = ("samtools view {} | grep -v 'AL:i:0' | awk ".format(tagBAM) +
               "'{{print $(NF-2)}}' | sed 's/BS:Z://g' | sort | uniq -c |"
               " awk 'BEGIN {{OFS=\"\t\"}};{{print $1,$2}}' > {}".format(idx))
        HTSUtils.PipedShellCall(Str)
    return idx


def GenBCIndexRef(idx, output="default"):
    pl("GenBCIndexRef. Input index: {}".format(idx))
    if(output == "default"):
        output = idx.split('.')[0] + '.ref.fasta'
    Str = "paste {0} {0} | sed 's:^:>:g' | tr '\t' '\n' > {1}".format(
          idx, output)
    HTSUtils.PipedShellCall(Str)
    return output


def GenerateFamilyHistochart(BCIdx, output="default"):
    if(output == "default"):
        output = '.'.join(BCIdx.split('.')[:-1]) + '.hist.txt'
    Str = ("cat {} | awk '{{print $1}}' | sort | uniq -c | ".format(BCIdx) +
           "awk 'BEGIN {{OFS=\"\t\"}};{{print $1,$2}}' | sort " +
           "-k1,1n > {}".format(output))
    pl("Command str: {}".format(Str))
    HTSUtils.PipedShellCall(Str)
    pl("Family size histochart: {}".format(output))
    return output


def getFamilySizeBAM(inputBAM, idx, output="default", passBC="default",
                     minFamSize=3):
    pl("getFamilySizeBAM. Input: {}".format(inputBAM))
    if(output == "default"):
        output = inputBAM.split('.')[0] + ".doubleFam.bam"
    if(passBC == "default"):
        passBC = inputBAM.split('.')[0] + ".doubleFam.lst"
    index = open(idx, "r")
    dictEntries = [line.split() for line in index]
    BarDict = {}
    for entry in dictEntries:
        BarDict[entry[1]] = entry[0]
    index.close()
    Reads = pysam.Samfile(inputBAM, "rb")
    outfile = pysam.Samfile(output, "wb", template=Reads)
    writeList = open(passBC, "w", 0)
    for record in Reads:
        Barcode = record.opt("BS")
        try:
            famSize = int(BarDict[Barcode])
        except KeyError:
            famSize = 0
        record.setTag("FM", famSize)
        if(int(famSize) >= minFamSize):
            writeList.write('{}\n'.format(Barcode))
            try:
                record.opt("BD")
            except KeyError:
                record.setTag("BD", 0)
        outfile.write(record)
    writeList.close()
    outfile.close()
    import uuid
    tempname = str(uuid.uuid4().get_hex().upper()[0:12]) + ".OMGZZZZ.tmp"
    string1 = "cat {0} | sort | uniq > {1} && mv {1} {0}".format(
        passBC, tempname)
    HTSUtils.PipedShellCall(string1)
    return output, passBC


def mergeBams(BAM1, BAM2, PT="default", outbam="default"):
    pl("mergeBams. BAM1: {}. BAM2: {}".format(BAM1, BAM2))
    if(PT == "default"):
        PT = "/mounts/bin/picard-tools"
    if(outbam == "default"):
        outbam = '.'.join(BAM1.split('.')[0:-1]) + '.merged.bam'
    command = ("java -jar {}/MergeSamFiles.jar I={}".format(PT, BAM1) +
               " I={} O={} SO=coordinate AS=true".format(BAM2, outbam))
    subprocess.check_call(shlex.split(command), shell=False)
    pl("Command string for merging was: {}".format(command))
    return outbam


def mergeBarcodes(reads1, reads2, outfile="default"):
    pl("mergeBarcodes. R1: {}. R2: {}".format(reads1, reads2))
    reader1 = pysam.Samfile(reads1, "rb")
    reader2 = pysam.Samfile(reads2, "rb")
    if(outfile == "default"):
        outfile = '.'.join(reads1.split('.')[0:-2]) + '.merged.bam'
    outSAM = pysam.Samfile(outfile, "wb", template=reader1)
    for entry1 in reader1:
        entry2 = reader2.next()
        assert entry1.qname == entry2.qname
        Barcode1 = entry1.opt("BS")
        Barcode2 = entry1.opt("BS")
        # print("Barcode 1: {}. Barcode 2: {}.".format(Barcode1,Barcode2))
        concatBarcode = Barcode1 + Barcode2
        # print("New barcode will be {}".format(concatBarcode))
        entry1.setTag("BS", concatBarcode)
        entry2.setTag("BS", concatBarcode)
        outSAM.write(entry1)
        outSAM.write(entry2)
    reader1.close()
    reader2.close()
    outSAM.close()
    return outfile


def pairedBarcodeTagging(
        fq1,
        fq2,
        bam,
        outBAMFile="default",
        suppBam="default"):
    if(outBAMFile == "default"):
        outBAMFile = '.'.join(bam.split('.')[0:-1]) + "tagged.bam"
    if(suppBam == "default"):
        suppBam = bam.split('.')[0] + '.2ndSupp.bam'
    pl("pairedBarcodeTagging. Fq: {}. outputBAM: {}".format(bam, outBAMFile))
    read1 = SeqIO.parse(fq1, "fastq")
    read2 = SeqIO.parse(fq2, "fastq")
    postFilterBAM = pysam.Samfile(bam, "rb")
    outBAM = pysam.Samfile(outBAMFile, "wb", template=postFilterBAM)
    suppBAM = pysam.Samfile(suppBam, "wb", template=postFilterBAM)
    for entry in postFilterBAM:
        if(entry.is_secondary or entry.flag > 2048):
            suppBAM.write(entry)
            continue
        if(not entry.is_paired):
            continue
        if(entry.is_read1):
            tempRead = read1.next()
            # print("Read desc: {}".format(tempRead.description))
        elif(entry.is_read2):
            tempRead = read2.next()
        descDict = BCFastq.GetDescriptionTagDict(tempRead.description)
        entry.setTag("FM", descDict["FM"])
        try:
            entry.setTag("BS", descDict["BS"])
        except KeyError:
            pl(("Dict: {}".format(descDict)))
            pl("Read: {}".format(entry))
            raise KeyError("Your fastq record is missing a BS tag.")
        try:
            if("Pass" in descDict["FP"]):
                entry.setTag("FP", 1)
            else:
                entry.setTag("FP", 0)
        except KeyError():
            pl(("Dict: {}".format(descDict)))
            pl("Read: {}".format(entry))
            raise KeyError("Your fastq record is missing an FP tag.")
        outBAM.write(entry)
    suppBAM.close()
    outBAM.close()
    postFilterBAM.close()
    return outBAMFile


def pairedFilterBam(inputBAM, passBAM="default",
                    failBAM="default", criteria="default",
                    deleteFailBam=False):
    """
    Filters out both reads in a pair based on a list of ","-separated criteria.
    Both reads must pass for the pair to be written
    Required: [sb]am file, name-sorted, keep unmapped reads,
    and remove supplementary and secondary alignments.
    criteria must be a string. If multiple filters are set,
    they must be comma-separated and have no spaces.
    """
    cStr = ("pairedFilterBam({}, passBAM='{}', ".format(inputBAM, passBAM) +
            "failBAM='{}', criteria='{}', delete".format(failBAM, criteria) +
            "FailBam='{}')".format(deleteFailBam))
    pl("Command required to reproduce this call: {}".format(cStr))
    if(criteria == "default"):
        raise NameError("Filter Failed: Criterion Not Set.")
    if(passBAM == "default"):
        passBAM = '.'.join(inputBAM.split('.')[0:-1])
        passBAM += ".{}P.bam".format(criteria)
    if(failBAM == "default"):
        failBAM = '.'.join(inputBAM.split('.')[0:-1])
        failBAM += ".{}F.bam".format(criteria)
    pl("pairedFilterBam. Input: {}. Pass: {}".format(inputBAM, passBAM))
    inBAM = pysam.Samfile(inputBAM, "rb")
    passFilter = pysam.Samfile(passBAM, "wb", template=inBAM)
    failFilter = pysam.Samfile(failBAM, "wb", template=inBAM)
    criteriaList = criteria.lower()
    pl("Criteria string is: {}".format(criteriaList))
    pl("Following criteria will be used: {}".format(
        str(criteriaList.split(','))))
    for read in inBAM:
        failed = False
        if(read.is_read1):
            read1 = read
            continue
        if(read.is_read2):
            read2 = read
            try:
                assert read1.qname == read2.qname  # Sanity check
            except AssertionError:
                pl("Failed sanity check. Is the alignment file name-sorted?")
                raise ThisIsMadness("Please inspect your bam file!")
            result = criteriaTest(read1, read2, filterSet=criteriaList)
            if(result is False):
                failed = True
                failFilter.write(read1)
                failFilter.write(read2)
                continue
            if(failed is True):
                continue
            passFilter.write(read1)
            passFilter.write(read2)
    passFilter.close()
    failFilter.close()
    inBAM.close()
    if(str(deleteFailBam).lower() == "true"):
        os.remove(failBAM)
    return passBAM, failBAM


def removeSecondary(inBAM, outBAM="default"):
    if(outBAM == "default"):
        outBAM = '.'.join(inBAM.split('.')[0:-1]) + '.2ndrm.bam'
    input = pysam.Samfile(inBAM, "rb")
    output = pysam.Samfile(outBAM, "wb", template=input)
    pl(("Attempting to remove secondary"))
    for entry in input:
        if(entry.is_secondary or entry.flag > 2048):
            continue
        else:
            output.write(entry)
    return outBAM


def Sam2Bam(insam, outbam):
    pl("Sam2Bam converting {} to {}".format(insam, outbam))
    output = open(outbam, 'w', 0)
    command_str = 'samtools view -Sbh {}'.format(insam)
    pl((command_str))
    subprocess.check_call(shlex.split(command_str), stdout=output, shell=False)
    return(command_str, outbam)


# Taken from Scrutils, by Jacob Durtschi
def SamtoolsBam2fq(bamPath, outFastqPath):
    pl("SamtoolsBam2fq converting {}".format(bamPath))
    # Build commands that will be piped
    samtoolsBam2fqCommand = ['samtools', 'bam2fq', bamPath]
    gzipCommand = ['gzip']

    # Open output fastq.gz file to be piped into
    outFastq = open(outFastqPath, 'w')

    # Call piped commands
    process1 = subprocess.Popen(samtoolsBam2fqCommand,
                                stdout=subprocess.PIPE, shell=False)
    process2 = subprocess.Popen(gzipCommand, stdin=process1.stdout,
                                stdout=outFastq, shell=False)
    # process1.stdout.close()
    # process1.wait()
    process2.wait()
    outFastq.flush()
    outFastq.close()
    # if processErr:

    return outFastqPath


def singleBarcodeTagging(fastq, bam, outputBAM="default", suppBam="default"):
    if(outputBAM == "default"):
        outputBAM = '.'.join(bam.split('.')[0:-1]) + "tagged.bam"
    if(suppBam == "default"):
        suppBam = bam.split('.')[0] + '.2ndSupp.bam'
    pl("singleBarcodeTagging. Fq: {}. outputBAM: {}".format(bam, outputBAM))
    reads = SeqIO.parse(fastq, "fastq")
    # inBAM = removeSecondary(args.bam_file) #Artefactual code
    postFilterBAM = pysam.Samfile(bam, "rb")
    suppBAM = pysam.Samfile(suppBam, "wb", template=postFilterBAM)
    outBAM = pysam.Samfile(outputBAM, "wb", template=postFilterBAM)
    for entry in postFilterBAM:
        if(entry.is_secondary or entry.flag > 2048):
            suppBAM.write(entry)
            continue
        else:
            try:
                tempRead = reads.next()
            except StopIteration:
                break
        descDict = BCFastq.GetDescriptionTagDict(tempRead.description)
        for key in descDict.keys():
            entry.setTag(key, descDict[key])
        if("Pass" not in descDict["FP"]):
            pl(("Standard filter for barcode failed! Val: {}".format(
                descDict["FP"])))
            pl(("Tags dictionary is {}".format(descDict)))
            raise ValueError("Something has gone wrong!")
        if("Pass" in descDict["FP"]):
            entry.tags = entry.tags + [("FP", 1)]
        else:
            entry.tags = entry.tags + [("FP", 0)]
        outBAM.write(entry)
    outBAM.close()
    postFilterBAM.close()
    return outputBAM


def SingleConsolidate(inbam, outbam="default", stringency=0.9):
    if(outbam == "default"):
        outbam = '.'.join(inbam.split('.')[0:-1]) + "consolidated.bam"
    pl("SingleConsolidate. Input: {}. Output: {}".format(inbam, outbam))
    inputHandle = pysam.Samfile(inbam, 'rb')
    outputHandle = pysam.Samfile(outbam, 'wb', template=inputHandle)
    workBC = ""
    Set = []
    for record in inputHandle:
        barcodeRecord = [
            tagSet for tagSet in record.tags if tagSet[0] == "BS"][0][1]
        # print("Working: {}. Current: {}.".format(workBC,barcodeRecord))
        # print("name of read with this barcode: {}".format(record.qname))
        # print("Working set: {}".format(Set))
        if(workBC == ""):
            workBC = barcodeRecord
            Set = []
            Set.append(record)
            continue
        elif(workBC == barcodeRecord):
            # print(record.qname)
            # print(workBC==barcodeRecord)
            # print("WorkBC: {}. Current: {}.".format(workBC,barcodeRecord))
            try:
                assert([tagSet for tagSet in Set[0].tags if tagSet[
                    0] == "BS"][0][1] == [
                    tagSet for tagSet in record.tags if tagSet[
                        0] == "BS"][0][1])
            except AssertionError:
                BS1 = [tagSet for tagSet in Set[
                    0].tags if tagSet[0] == "BS"][0][1]
                BS2 = [tagSet for tagSet in record.tags if tagSet[
                    0] == "BS"][0][1]
                pl("BS1: {}. BS2: {}".format(BS1, BS2))
                raise AssertionError("Well, there you go.")
            Set.append(record)
            continue
        elif(workBC != barcodeRecord):
            mergeRec, success = compareRecs(Set, stringency=stringency)
            if(success is True):
                outputHandle.write(mergeRec)
            Set = []
            Set.append(record)
            workBC = barcodeRecord
            continue
        else:
            raise RuntimeError("This code should be unreachable")
    inputHandle.close()
    outputHandle.close()
    return outbam


def singleCriteriaTest(read, filter="default"):
    list = "adapter complexity editdistance family ismapped qc".split(' ')

    if(filter == "default"):
        pl(("List of valid filters: {}".format(', '.join(list))))
        raise ValueError("Filter must be set!.")

    if(filter not in list):
        raise ValueError("Filter not supported. The list is: {}".format(list))

    if(filter == "adapter"):
        ALloc1 = [i for i, j in enumerate(read.tags) if j[0] == "FP"][0]

        try:
            ALValue1 = int(read.tags[ALloc1][1])
        except IndexError:
            ALValue1 = 0
        if(ALValue1 == 0):
            return False

    if(filter == "complexity"):
        a = "AAAAAAAAAA"
        t = "TTTTTTTTTT"
        c = "CCCCCCCCCC"
        g = "GGGGGGGGGG"
        rt = str(read.tags)
        if(a in rt or t in rt or g in rt or c in rt):
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


def singleFilterBam(inputBAM, passBAM="default",
                    failBAM="default", criteria="default"):
    """
    Filters out both reads in a pair based on a list of
    comma-separated criteria. Both reads must pass to be written.
    Required: [SB]AM file, coordinate-sorted,
    supplementary and secondary alignments removed,unmapped reads retained.
    """
    if(criteria == "default"):
        raise NameError("Filter Failed: Criterion Not Set.")
    if(passBAM == "default"):
        passBAM = '.'.join(inputBAM.split('.')[0:-1])
        passBAM += ".{}P.bam".format(criteria)
    if(failBAM == "default"):
        failBAM = '.'.join(inputBAM.split('.')[0:-1])
        failBAM += ".{}F.bam".format(criteria)
    pl("singleFilterBam. Input: {}. Pass: {}".format(inputBAM, passBAM))
    inBAM = pysam.Samfile(inputBAM, "rb")
    passFilter = pysam.Samfile(passBAM, "wb", template=inBAM)
    failFilter = pysam.Samfile(failBAM, "wb", template=inBAM)
    criteriaList = criteria.lower().split(',')
    for i, entry in enumerate(criteriaList):
        pl(("Criteria #{} is \"{}\"".format(i, entry)))
    for read in inBAM:
        failed = False
        for criterion in criteriaList:
            result = singleCriteriaTest(read, filter=criterion)
            if(result is False):
                failFilter.write(read)
                failed = True
                break
            else:
                continue
        if(failed is True):
            continue
        passFilter.write(read)
    passFilter.close()
    failFilter.close()
    inBAM.close()
    return passBAM, failBAM


def splitBAMByReads(BAM, read1BAM="default", read2BAM="default"):
    pl("splitBAMByReads. Input: {}".format(BAM))
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
