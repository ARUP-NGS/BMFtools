#cython: boundscheck=False
import logging
import os
import shlex
import subprocess
import copy

from Bio import SeqIO
import Bio
from Bio.SeqRecord import SeqRecord
import cython
cimport cython
import numpy as np
cimport numpy as np

from utilBMF.HTSUtils import printlog as pl, ThisIsMadness
from utilBMF.HTSUtils import PipedShellCall
from utilBMF import HTSUtils

letterNumDict = {}
letterNumDict['A'] = 0
letterNumDict['C'] = 1
letterNumDict['G'] = 2
letterNumDict['T'] = 3
letterNumDict[0] = 'A'
letterNumDict[1] = 'C'
letterNumDict[2] = 'G'
letterNumDict[3] = 'T'

"""
Contains various utilities for working with barcoded fastq files.
TODO: Convert fastq reading to pysam.FastqFile (reading)
and to writing as a string
FastqStr = "\n".join(["@" + c.name + c.comment, c.sequence , "+", c.quality])
FastqSetStr = "\n".join(["\n".join(["@" + c.name + c.comment, c.sequence ,
    "+", c.quality]) for c in FastqRecordsList])
GzipFastqSetStr = zlib.compress(FastqSetStr)
This will add support for gzipped FastqFiles as well as some
significant performance increases.
"""


def BarcodeSort(inFastq, outFastq="default"):
    pl("Sorting {} by barcode sequence.".format(inFastq))
    if(outFastq == "default"):
        outFastq = '.'.join(inFastq.split('.')[0:-1]) + '.BS.fastq'
    BSstring = ("cat " + inFastq + " | paste - - - - | grep 'Pass' | sed "
                "'s: #G~:\t#G~:g' |  awk 'BEGIN {{FS=OFS=\"\t\"}};{{print"
                " $3,$0}}' | sort -k1,1 | cut -f2- | "
                "sed 's:\t#G~: #G~:g' | tr '\t' '\n' > " + outFastq)
    PipedShellCall(BSstring)
    # pl("Command: {}".format(BSstring.replace(
    #    "\t", "\\t").replace("\n", "\\n")))
    return outFastq


@cython.locals(stringency=cython.float, hybrid=cython.bint,
               famLimit=cython.int, keepFails=cython.bint,
               Success=cython.bint, PASS=cython.bint, frac=cython.float)
def compareFastqRecords(R, stringency=0.9, hybrid=False, famLimit=100,
                        keepFails=True):
    """
    Compares the fastq records to create a consensus sequence (if it
    passes a filter)
    """
    try:
        famLimit = int(famLimit)
    except ValueError:
        pl("famLimit arg must be integer. Set to default: 100.")
    if(len(R) > famLimit):
        logging.debug(
            "Read family - {} with {} members was capped at {}. ".format(
                R[0], len(R), famLimit))
        R = R[:famLimit]
        # Add debugging logging for this step.
    seqs = [str(record.seq) for record in R]
    maxScore = 0
    Success = False
    for seq in seqs:
        # print("Seq: {}".format(str(seq)))
        numEq = sum(str(seq) == str(seqItem) for seqItem in seqs)
        if(numEq > maxScore):
            maxScore = numEq
            finalSeq = str(seq)
    frac = numEq * 1.0 / len(R)
    PASS = frac > stringency
    # print("Fraction {}. Stringency: {}. Pass? {}.".format(
    # frac,stringency,PASS))
    if(PASS):
        Success = True
    elif(frac < 0.5):
        Success = False
    elif(hybrid is True):
        return compareFastqRecordsInexactNumpy(R)
    consolidatedRecord = SeqRecord(seq=finalSeq, id=R[0].id,
                                   letter_annotations=R[0].letter_annotations,
                                   name=R[0].name,
                                   description=R[0].description)
    if("Fail" in GetDescTagValue(consolidatedRecord.description, "FP")):
        Success = False
    probs = np.multiply(len(R), R[0].letter_annotations['phred_quality'])
    if(np.any(np.less(probs, 1))):
        Success = False
    if(np.any(np.greater(probs, 93))):
        consolidatedRecord.description += " #G~PV=" + ",".join(
            probs.astype(str))
    probs[probs <= 0] = 93
    probs[probs > 93] = 93
    consolidatedRecord.letter_annotations[
        'phred_quality'] = probs
    return consolidatedRecord, Success


@cython.locals(Success=cython.bint)
def compareFastqRecordsInexactNumpy(R):
    """
    Calculates the most likely nucleotide
    at each position and returns the joined record.
    """

    def dAccess(x):
        return letterNumDict[x]
    dAccess = np.vectorize(dAccess)
    seqs = [str(record.seq) for record in R]
    stackArrays = tuple([np.char.array(s, itemsize=1) for s in seqs])
    seqArray = np.vstack(stackArrays)
    # print(repr(seqArray))
    cdef np.ndarray[np.int, ndim = 1] quals = np.array([
        record.letter_annotations['phred_quality'] for record in R])
    cdef np.ndarray[np.int, ndim = 1] qualA = copy.copy(quals)
    cdef np.ndarray[np.int, ndim = 1] qualC = copy.copy(quals)
    cdef np.ndarray[np.int, ndim = 1] qualG = copy.copy(quals)
    cdef np.ndarray[np.int, ndim = 1] qualT = copy.copy(quals)
    qualA[seqArray != "A"] = 0
    qualA = np.sum(qualA, 0)
    qualC[seqArray != "C"] = 0
    qualC = np.sum(qualC, 0)
    qualG[seqArray != "G"] = 0
    qualG = np.sum(qualG, 0)
    qualT[seqArray != "T"] = 0
    qualT = np.sum(qualT, 0)
    cdef np.ndarray[np.int, ndim = 2] qualAllSum = np.vstack(
        [qualA, qualC, qualG, qualT])
    Success = True
    newSeq = "".join(
        np.apply_along_axis(dAccess, 0, np.argmax(qualAllSum, 0)))
    cdef np.ndarray[np.int, ndim = 1] MaxPhredSum = np.amax(
        qualAllSum, 0)  # Avoid calculating twice.
    cdef np.ndarray[np.int, ndim = 1] phredQuals = np.subtract(
        np.multiply(2, MaxPhredSum), np.sum(qualAllSum, 0))
    phredQuals[phredQuals == 0] = 93
    phredQuals[phredQuals < 0] = 0
    consolidatedRecord = SeqRecord(
        seq=newSeq,
        id=R[0].id,
        name=R[0].name,
        description=R[0].description)
    if(np.any(np.greater(phredQuals, 93))):
        consolidatedRecord.description += (" #G~PV=" +
                                           ",".join(phredQuals.astype(str)))
    phredQuals[phredQuals > 93] = 93
    consolidatedRecord.letter_annotations[
        'phred_quality'] = phredQuals.tolist()
    if("Fail" in GetDescTagValue(consolidatedRecord.description, "FP")):
        Success = False
    # Checking for a strange failure on the part of the package
    # to correctly assign the phred score.
    try:
        RealQuals = [int(i) for i in GetDescTagValue(
            consolidatedRecord.description, "PV").split(',')]
    except KeyError:
        return consolidatedRecord, Success
    try:
        assert np.all([q == r or (r > q and q == 93) for q, r in zip(
            consolidatedRecord.letter_annotations['phred_quality'],
            RealQuals)])
    except AssertionError:
        pl("phred quality seemed to have been jumbled just a bit - fixing!")
        newQuals = []
        for q, r in zip(consolidatedRecord.letter_annotations['phred_quality'],
                        RealQuals):
            if(q > r):
                newQuals.append(q)
            else:
                newQuals.append(r)
            newQuals = [i if i <= 93 else 93 for i in newQuals]
        try:
            consolidatedRecord.letter_annotations['phred_quality'] = newQuals
        except TypeError:
            pl(repr(newQuals) + " and length {}".format(len(newQuals)))
            raise ThisIsMadness("Problem with quality strings.")
    return consolidatedRecord, Success


@cython.locals(gzip=cython.bint, bLen=cython.int)
def FastqPairedShading(fq1,
                       fq2,
                       indexfq="default",
                       outfq1="default",
                       outfq2="default",
                       gzip=True):
    """
    Tags fastqs with barcodes from an index fastq.
    TODO: change the write-out to concatenate strings before writing.
    """
    pl("Now beginning fastq marking: Pass/Fail and Barcode")
    if(indexfq == "default"):
        raise ValueError("For an i5/i7 index ")
    if(outfq1 == "default"):
        outfq1 = ('.'.join(
            fq1.split('.')[0:-1]) + '.shaded.fastq').split('/')[-1]
    if(outfq2 == "default"):
        outfq2 = ('.'.join(
            fq2.split('.')[0:-1]) + '.shaded.fastq').split('/')[-1]
    pl("Output fastqs: {}, {}.".format(outfq1, outfq2))
    inFq1 = SeqIO.parse(fq1, "fastq")
    inFq2 = SeqIO.parse(fq2, "fastq")
    outFqHandle1 = open(outfq1, "w")
    outFqHandle2 = open(outfq2, "w")
    inIndex = SeqIO.parse(indexfq, "fastq")
    for read1 in inFq1:
        read2 = inFq2.next()
        indexRead = inIndex.next()
        tempBar, bLen = indexRead.seq, len(indexRead.seq) * 5 / 6
        # bLen - 10 of 12 in a row, or 5/6. See Loeb, et al.
        # This is for removing low complexity reads
        # print("bLen is {}".format(bLen))
        if(("N" in tempBar or "A"*bLen in tempBar
                or "C"*bLen in tempBar or "G"*bLen in tempBar
                or "T"*bLen in tempBar)):
            pl("Failing barcode for read {} is {} ".format(indexRead,
                                                           tempBar),
               level=logging.DEBUG)
            read1.description += " #G~FP=IndexFail #G~BS=" + str(indexRead.seq)
            read2.description += " #G~FP=IndexFail #G~BS=" + str(indexRead.seq)
            SeqIO.write(read1, outFqHandle1, "fastq")
            SeqIO.write(read2, outFqHandle2, "fastq")
        else:
            read1.description += " #G~FP=IndexPass #G~BS=" + str(indexRead.seq)
            read2.description += " #G~FP=IndexPass #G~BS=" + str(indexRead.seq)
            SeqIO.write(read1, outFqHandle1, "fastq")
            SeqIO.write(read2, outFqHandle2, "fastq")
    outFqHandle1.close()
    outFqHandle2.close()
    if(gzip is True):
        from subprocess import check_call
        check_call(['gzip', fq1], shell=False)
        check_call(['gzip', fq2], shell=False)
    return outfq1, outfq2


def FastqSingleShading(fq,
                       indexfq="default",
                       outfq="default",
                       gzip=True):
    pl("Now beginning fastq marking: Pass/Fail and Barcode")
    if(indexfq == "default"):
        raise ValueError("For an i5/i7 index ")
    if(outfq == "default"):
        outfq = '.'.join(fq.split('.')[0:-1]) + '.shaded.fastq'
    inFq1 = SeqIO.parse(fq, "fastq")
    outFqHandle1 = open(outfq, "w")
    inIndex = SeqIO.parse(indexfq, "fastq")
    for read1 in inFq1:
        indexRead = inIndex.next()
        if("N" in indexRead.seq):
            read1.description += " #G~FP=IndexFail #G~BS=" + indexRead.seq
            SeqIO.write(read1, outFqHandle1, "fastq")
        else:
            read1.description += " #G~FP=IndexPass #G~BS=" + indexRead.seq
            SeqIO.write(read1, outFqHandle1, "fastq")
    outFqHandle1.close()
    if(gzip is True):
        from subprocess import check_call
        check_call(['gzip', fq], shell=False)
    return


def HomingSeqLoc(fq, homing, bar_len=12):
    pl("Now beginning HomingSeqLoc.")
    InFastq = SeqIO.parse(fq, "fastq")
    Tpref = '.'.join(fq.split('.')[0:-1])
    Prefix = Tpref.split('/')[-1]
    StdFilename = Prefix + '.{}.fastq'.format("homing" + str(bar_len))
    ElseFilename = Prefix + '.else.fastq'
    ElseLoc = Prefix + '.else.supp'
    StdFastq = open(StdFilename, 'w', 0)  # Homing at expected Location
    ElseFastq = open(ElseFilename, 'w', 0)
    ElseLocations = open(ElseLoc, 'w', 0)
    for read in InFastq:
        seq = str(read.seq)
        if(seq.find(homing) == -1):
            read.description += " #G~FP=HomingFail"
            SeqIO.write(read, StdFastq, "fastq")
        elif(seq[bar_len:bar_len + len(homing)] == homing):
            read.description += " #G~FP=HomingPass"
            SeqIO.write(read, StdFastq, "fastq")
        else:
            read.description = " #G~FP=HomingFail"
            SeqIO.write(read, StdFastq, "fastq")
            ElseLocations.write(repr(seq.find(homing)) + "\t" +
                                read.name + "\n")
    StdFastq.close()
    ElseFastq.close()
    ElseLocations.close()
    return StdFilename, ElseFilename


def fastq_sort(in_fastq, out_fastq):
    pl("Now beginning fastq_sort.")
    import subprocess
    outfile = open(out_fastq, 'w')
    command_str = ('cat ' + in_fastq + ' | paste - - - - | '
                   'sort -k1,1 -t " " | tr "\t" "\n"')
    subprocess.check_call(command_str, stdout=open(outfile, "w"), shell=True)
    outfile.close()
    return(command_str)


def FastqRegex(fq, string, matchFile="default", missFile="default"):
    from subprocess import check_call
    if(matchFile == "default"):
        matchFile = ('.'.join(fq.split(
                     '.')[0:-1]) + '.match.fastq').split('/')[-1]
    if(missFile == "default"):
        missFile = (
            '.'.join(fq.split('.')[0:-1]) + '.miss.fastq').split('/')[-1]
    CommandStr = ("cat " + fq + " | paste - - - - | grep '" +
                  string + "' | tr '\t' '\n' > " + matchFile)
    check_call(CommandStr, shell=True)
    CommandStr2 = ("cat {} | paste - - - - | grep ".format(fq) +
                   "-v '{}' | tr '\t' '\n' > {}".format(
        string,
        missFile))
    check_call(CommandStr2, shell=True)
    return(CommandStr, CommandStr2, matchFile, missFile)


def fastx_trim(infq, outfq, n):
    pl("Now beginning fastx_trimmer.")
    command_str = ['fastx_trimmer', '-l', str(n), '-i', infq, '-o', outfq]
    pl(command_str)
    subprocess.check_call(command_str)
    return(command_str)


def GetDescTagValue(readDesc, tag="default"):
    """
    Gets the value associated with a given tag.
    If a SeqRecord object rather than a string
    of read description (string), then the seq
    attribute is used instead of the description.
    """
    if(tag == "default"):
        raise ValueError("A tag must be specified!")
    try:
        return GetDescriptionTagDict(readDesc)[tag]
    except KeyError:
        # pl("Tag {} is not available in the description.".format(tag))
        # pl("Description: {}".format(readDesc))
        raise KeyError("Invalid tag.")


def GetDescriptionTagDict(readDesc):
    """Returns a set of key/value pairs in a dictionary for """
    tagSetEntries = [i.strip().split("=") for i in readDesc.split("#G~")][1:]
    tagDict = {}
    try:
        for pair in tagSetEntries:
            tagDict[pair[0]] = pair[1]
    except IndexError:
        pl("A value is stored with the #G~ tag which doesn't contain an =.")
        pl("tagSetEntries: {}".format(tagSetEntries))
        raise IndexError("Check that fastq description meets specifications.")
    # pl("Repr of tagDict is {}".format(tagDict))
    except TypeError:
        pl("tagSetEntries: {}".format(tagSetEntries))
        ThisIsMadness("YOU HAVE NO CHANCE TO SURVIVE MAKE YOUR TIME")
    return tagDict


def GenerateShadesIndex(indexFastq, index_file="default"):
    from subprocess import check_call
    if(index_file == "default"):
        index_file = '.'.join(indexFastq.split('.')[0:-1]) + ".barIdx"
    commandStr = ("cat {} | paste - - - - | cut -f2 | ".format(indexFastq) +
                  "sort | uniq -c | awk 'BEGIN {{OFS=\"\t\"}};{{print $1"
                  ",$2}}' > {}".format(index_file))
    check_call(commandStr, shell=True)
    return index_file


# Replaces GenerateOnePairFastqBarcodeIndex
def GenerateSingleBarcodeIndex(tags_file, index_file="default"):
    import collections
    pl("Now beginning GenerateSingleBarcodeIndex for {}.".format(tags_file))
    if(index_file == "default"):
        index_file = '.'.join(tags_file.split('.')[0:-1]) + ".barIdx"
    index_handle = open(index_file, "w")
    inFq = SeqIO.parse(tags_file, "fastq")
    barcodeList = [GetDescTagValue(read.description, "BS") for read in inFq]
    counts = collections.Counter(barcodeList)
    for key in counts.keys():
        index_handle.write("{}\t{}\n".format(counts[key], key))
    return index_file


@cython.locals(deleteInFqs=cython.bint, minFamSize=cython.int)
def GetFamilySizePaired(
        trimfq1,
        trimfq2,
        BarcodeIndex,
        outfq1="default",
        outfq2="default",
        singlefq1="default",
        singlefq2="default",
        minFamSize=1,
        deleteInFqs=True):
    pl("Running GetFamilySizeSingle for {}, {}.".format(trimfq1, trimfq2))
    infq1 = SeqIO.parse(trimfq1, "fastq")
    infq2 = SeqIO.parse(trimfq2, "fastq")
    if(outfq1 == "default"):
        outfq1 = '.'.join(trimfq1.split('.')[0:-1] + ['fam', 'fastq'])
    if(outfq2 == "default"):
        outfq1 = '.'.join(trimfq2.split('.')[0:-1] + ['fam', 'fastq'])
    outfqBuffer1 = open(outfq1, "w")
    outfqBuffer2 = open(outfq2, "w")
    if(singlefq1 == "default"):
        singlefq1 = '.'.join(
            trimfq1.split('.')[
                0:-1]) + ".lonely.hearts.club.band.fastq"
    if(singlefq2 == "default"):
        singlefq2 = '.'.join(
            trimfq2.split('.')[
                0:-1]) + ".lonely.hearts.club.band.fastq"
    singlefqBuffer1 = open(singlefq1, "w")
    singlefqBuffer2 = open(singlefq2, "w")
    TotalReads, ReadsWithFamilies = 0, 0
    index = open(BarcodeIndex, "r")
    dictEntries = [line.split() for line in index]
    BarDict = {}
    for entry in dictEntries:
        BarDict[entry[1]] = entry[0]
    for read in infq1:
        read2 = infq2.next()
        # index.seek(0)
        TotalReads += 1
        readTag1 = GetDescTagValue(read.description, "BS")
        newRead1 = read
        newRead2 = read2
        try:
            famSize = BarDict[readTag1]
        except KeyError:
            famSize = "0"
        newRead1.description = read.description + " #G~FM=" + famSize
        newRead2.description = read2.description + " #G~FM=" + famSize
        # print("famSize = _{}_".format(str(famSize)))
        # print("The value of this comparison to 1 is {}".format(
        # str(famSize=="1")))
        if(int(famSize) <= int(minFamSize)):
            #  print("Hey, I found a singleton")
            SeqIO.write(newRead1, singlefqBuffer1, "fastq")
            SeqIO.write(newRead2, singlefqBuffer2, "fastq")
        else:
            #  print("Hey, I found a read with a family!")
            ReadsWithFamilies += 1
            SeqIO.write(newRead1, outfqBuffer1, "fastq")
            SeqIO.write(newRead2, outfqBuffer2, "fastq")
    outfqBuffer1.close()
    outfqBuffer2.close()
    singlefqBuffer1.close()
    singlefqBuffer2.close()
    if(deleteInFqs is True):
        os.remove(trimfq1)
        os.remove(trimfq2)
    return ([outfq1, outfq2],
            [singlefq1, singlefq2],
            TotalReads,
            ReadsWithFamilies)


def GetFamilySizeSingle(
        trimfq,
        BarcodeIndex,
        outfq="default",
        singlefq="default",
        minFamSize=1):
    pl("Running GetFamilySizeSingle for {}.".format(trimfq))
    infq = SeqIO.parse(trimfq, "fastq")
    if(outfq == "default"):
        outfq = '.'.join(trimfq.split('.')[0:-1]) + ".fam.fastq"
    outfqBuffers = open(outfq, "w", 0)
    index = open(BarcodeIndex, "r")
    if(singlefq == "default"):
        singlefq = '.'.join(
            trimfq.split('.')[
                0:-1]) + ".lonely.hearts.club.band.fastq"
    singlefqBuffer = open(singlefq, "w", 0)
    TotalReads, ReadsWithFamilies = 0, 0
    dictEntries = [line.split() for line in index]
    BarDict = {}
    for entry in dictEntries:
        BarDict[entry[1]] = entry[0]
    for read in infq:
        # index.seek(0)
        TotalReads += 1
        # print("description is {}".format(read.description))
        readTag = GetDescTagValue(read.description, "BS")
        newRead = read
        # print("readTag is _" + readTag + "_")
        try:
            famSize = BarDict[readTag]
        except KeyError:
            famSize = 0
        newRead.description = read.description + " #G~FM=" + str(famSize)
        if(str(famSize) == "0"):
            continue
        if(int(famSize) < minFamSize):
            #  print("Hey, I found a singleton")
            SeqIO.write(newRead, singlefqBuffer, "fastq")
        else:
            #  print("Hey, I found a read with a family!")
            ReadsWithFamilies += 1
            SeqIO.write(newRead, outfqBuffers, "fastq")
    outfqBuffers.close()
    singlefqBuffer.close()
    return outfq, TotalReads, ReadsWithFamilies


def LighterCallPaired(fq1, fq2, kmer="default",
                      captureSize="default",
                      alpha="default"):
    """
    Calls Lighter on both sets of fastq files. Default kmer of 20.
    Capture size is required. If alpha is unset, lighter will infer the value.
    """
    if(kmer == "default"):
        kmer = 20
        pl("kmer not set - default of 20 used.")
    outfq1 = ".".join(fq1.split(".")[0:-1] + ["cor", "fq"])
    outfq2 = ".".join(fq2.split(".")[0:-1] + ["cor", "fq"])
    if(isinstance(captureSize, str)):
        raise ThisIsMadness("Capture size must be set for error correction.")
    if(isinstance(alpha, str)):
        pl("Alpha constant not provided. Lighter will infer it from captureS"
           "ize, the number of reads, and the length of each.")
        commandArray = shlex.split("lighter -r {} -r {} -K {} {}".format(
            fq1, fq2, kmer, captureSize))
        subprocess.check_call(commandArray, shell=False)
    else:
        commandArray = shlex.split("lighter -r {} -r {} -K {} {} {}".format(
            fq1, fq2, kmer, captureSize, alpha))
    return outfq1, outfq2


def mergeBarcodes(fq1, fq2, out1="default", out2="default"):
    if(out1 == "default"):
        out1 = '.'.join(fq1.split('.')[0:-1]) + '.mergeBC.fastq'
    if(out2 == "default"):
        out2 = '.'.join(fq2.split('.')[0:-1]) + '.mergeBC.fastq'
    reads1 = SeqIO.parse(fq1, "fastq")
    reads2 = SeqIO.parse(fq2, "fastq")
    outFq1 = open(out1, "w")
    outFq2 = open(out2, "w")
    while True:
        try:
            read1 = reads1.next()
            read2 = reads2.next()
            read1Vals = GetDescriptionTagDict(read1.description)
            read2Vals = GetDescriptionTagDict(read2.description)
            if("Fail" in read1Vals["FP"] + read2Vals["FP"]):
                continue
            newBarcode = read1Vals["BS"] + read2Vals["BS"]
            read1.description = read1.description.replace(
                "BS=" + read1Vals["BS"], "BS=" + newBarcode)
            read2.description = read2.description.replace(
                "BS=" + read2Vals["BS"], "BS=" + newBarcode)
            SeqIO.write(read1, outFq1, "fastq")
            SeqIO.write(read2, outFq2, "fastq")
        except StopIteration:
            break
    reads1.close()
    reads2.close()
    return out1, out2


def mergeSequencesFastq(fq1, fq2, output="default"):
    pl("mergeSequencesFastq for {} and {}".format(fq1, fq2))
    if(output == "default"):
        output = fq1.split('.')[0] + '.merged.fastq'
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
    reads1 = SeqIO.parse(fq1, "fastq")
    reads2 = SeqIO.parse(fq2, "fastq")
    outFastq = open(output, "w", 0)
    while True:
        try:
            read1 = reads1.next()
            read2 = reads2.next()
            outread = read1
            outread = SeqRecord(
                Seq(str(read1.seq) + str(read2.seq), "fastq"),
                id=read1.id, description=read1.description)
            outread.letter_annotations[
                'phred_quality'] = read1.letter_annotations[
                'phred_quality'] + read2.letter_annotations['phred_quality']
            SeqIO.write(outread, outFastq, "fastq")
        except StopIteration:
            break
    reads1.close()
    reads2.close()
    outFastq.close()
    return output


def PairFastqBarcodeIndex(taggedFile1, taggedFile2, index_file="default"):
    pl("Now beginning GenerateFullFastqBarcodeIndex for {} and {}.".format(
        taggedFile1, taggedFile2))
    from subprocess import check_call
    if(index_file == "default"):
        index_file = '.'.join(taggedFile1.split('.')[0:-1]) + ".barIdx"
    cmd = "cat {} {} | sed 's: #G~BS=:\t:g' | paste - - - - | awk ".format(
        taggedFile1, taggedFile2)
    cmd += "'BEGIN {{FS=\"\t\"}};{{print $3}}' | sort | uniq -c | awk 'BEGIN "
    cmd += "{{OFS=\"\t\"}};{{print $1,$2}}' | sort -k1,1n > {}".format(
        index_file)
    pl("CommandStr = {}".format(cmd.replace("\t", "\\t")))
    check_call(cmd, shell=True)
    return index_file


@cython.locals(stringency=cython.float,
               numpy=cython.bint, keepFailedPairs=cython.bint)
def pairedFastqConsolidate(fq1, fq2, outFqPair1="default",
                           outFqPair2="default",
                           stringency=0.9,
                           numpy=False,
                           keepFailedPairs=True):
    pl("Now running pairedFastqConsolidate on {} and {}.".format(fq1, fq2))
    pl(("Command required to duplicate this action:"
        " pairedFastqConsolidate('{}', '{}', outFqPair1".format(fq1, fq2)) +
        "='{}', outFqPair2='{}', ".format(outFqPair1, outFqPair2) +
        "stringency='{}', numpy='{}')".format(stringency, numpy))
    if(outFqPair1 == "default"):
        outFqPair1 = '.'.join(fq1.split('.')[0:-1]) + 'cons.fastq'
    if(outFqPair2 == "default"):
        outFqPair2 = '.'.join(fq2.split('.')[0:-1]) + 'cons.fastq'
    inFq1 = SeqIO.parse(fq1, 'fastq')
    inFq2 = SeqIO.parse(fq2, 'fastq')
    outputHandle1 = open(outFqPair1, 'w')
    outputHandle2 = open(outFqPair2, 'w')
    workingBarcode = ""
    workingSet1 = []
    workingSet2 = []
    for fqRec in inFq1:
        try:
            fqRec2 = inFq2.next()
        except ValueError:
            print("Trying to read in from " + fq2)
            print(fqRec)
            print inFq2.next()
            raise ValueError("I am confused.")
        bc4fq1 = GetDescTagValue(fqRec.description, "BS")
        # Originally removing reads with family size <2, since one pair could
        # have more than the other, it's important that I keep these reads in
        # and filter them from the BAM filea
        if(workingBarcode == ""):
            workingBarcode = bc4fq1
            workingSet1 = []
            workingSet1.append(fqRec)
            workingSet2 = []
            workingSet2.append(fqRec2)
            continue
        elif(workingBarcode == bc4fq1):
            workingSet1.append(fqRec)
            workingSet2.append(fqRec2)
            continue
        elif(workingBarcode != bc4fq1):
            # pl("About to collapse a family.")
            if(numpy is False):
                mergedRecord1, success1 = compareFastqRecords(
                    workingSet1, stringency=float(stringency))
                mergedRecord2, success2 = compareFastqRecords(
                    workingSet2, stringency=float(stringency))
            else:
                mergedRecord1, success1 = compareFastqRecords(
                    workingSet1, hybrid=True)
                mergedRecord2, success2 = compareFastqRecords(
                    workingSet2, hybrid=True)
            if(success1 is False):
                mergedRecord1.description.replace("Pass", "Fail")
            if(success2 is False):
                mergedRecord2.description.replace("Pass", "Fail")
            # Remove the " bad_prefix" added to the fastq name by Lighter.
            # Compatibility purposes.
            if(keepFailedPairs is True):
                SeqIO.write(mergedRecord2, outputHandle2, "fastq")
                SeqIO.write(mergedRecord1, outputHandle1, "fastq")
            elif("Fail" not in (mergedRecord1.description +
                                mergedRecord2.description)):
                SeqIO.write(mergedRecord2, outputHandle2, "fastq")
                SeqIO.write(mergedRecord1, outputHandle1, "fastq")
            workingSet1 = []
            workingSet1.append(fqRec)
            workingSet2 = []
            workingSet2.append(fqRec2)
            workingBarcode = bc4fq1
            continue
        else:
            raise RuntimeError(
                "No idea what's going on. This code should be unreachable")
    inFq1.close()
    inFq2.close()
    outputHandle1.close()
    outputHandle2.close()
    return outFqPair1, outFqPair2


def renameReads(fq1, fq2, outfq1="default", outfq2="default"):
    """
    Requires barcode-sorted, consolidated fastq files,
    filtered to only the shared barcode families
    """
    if(outfq1 == "default"):
        outfq1 = fq1.split('.')[0] + '.cons.R1fastq'
    if(outfq2 == "default"):
        outfq2 = fq2.split('.')[0] + '.cons.R2fastq'
    infq1 = SeqIO.parse(fq1, "fastq")
    infq2 = SeqIO.parse(fq2, "fastq")
    outfqhandle1 = open(outfq1, "w")
    outfqhandle2 = open(outfq2, "w")
    for read1 in infq1:
        read2 = infq2.next()
        read2.id = read1.id
        read2.description = ' '.join(
            read1.description.split()[1:]).replace('1:N', '2:N')
        SeqIO.write(read1, outfqhandle1, "fastq")
        SeqIO.write(read2, outfqhandle2, "fastq")
    outfqhandle1.close()
    outfqhandle2.close()
    infq1.close()
    infq2.close()
    return outfq1, outfq2


def ShadesRescuePaired(singlefq1,
                       singlefq2,
                       toBeRescuedFq1="default",
                       toBeRescuedFq2="default",
                       fqToAlign="default",
                       appendFq1="default",
                       appendFq2="default",
                       index="default",
                       minFamSize=3,
                       maxMismatches=2):
    """
    This extracts a list of all barcodes which have families
    from the index file provided.
    It then creates a custom reference for those barcodes.
    It then creates a dummy fastq which it then aligns to the custom
    reference.
    Finally, these reads have their barcodes changed to match the family's,
    and a BD:i: tag is added to specify how far removed it was.
    This function has also not been fully written
    """
    if(appendFq1 == "default" or appendFq2 == "default"):
        raise ValueError("Fastq for appending rescued reads required.")
    if(index == "default"):
        raise ValueError("Index for initial family sizes required.")
    ind = open(index, "r")
    dictEntries = [line.split() for line in ind]
    BarDict = {}
    for entry in dictEntries:
        if(int(entry[0]) >= minFamSize):
            BarDict[entry[1]] = entry[0]
    readFq1 = SeqIO.parse(singlefq1, "fastq")
    readFq2 = SeqIO.parse(singlefq2, "fastq")
    readToRescue1 = open(toBeRescuedFq1, "w")
    readToRescue2 = open(toBeRescuedFq2, "w")
    for read1 in readFq1:
        read2 = readFq2.next()
        try:
            GetDescTagValue(read1.desc, "FM")
        except KeyError:
            SeqIO.write(read1, readToRescue1, "fastq")
            SeqIO.write(read2, readToRescue2, "fastq")
    raise ThisIsMadness("This rescue function has not been written.")


def singleFastqConsolidate(fq, outFq="default", stringency=0.9):
    if(outFq == "default"):
        outFq = '.'.join(fq.split('.')[0:-1]) + 'cons.fastq'
    inFq = SeqIO.parse(fq, 'fastq')
    outputHandle = open(outFq, 'w')
    workingBarcode = ""
    workingSet = []
    for fqRec in inFq:
        bc4fq = GetDescTagValue(fqRec.description, "BS")
        # print("Working barcode: {}. Current barcode: {}.".format(
        #   workingBarcode,barcodeRecord))
        # print("name of read with this barcode: {}".format(record.qname))
        # print("Working set: {}".format(workingSet))
        bLen = len(bc4fq) * 5 / 6
        # bLen - 10 of 12 in a row, or 5/6. See Loeb, et al.
        # This is for removing low complexity reads
        if(("N" in bc4fq or "A"*bLen in bc4fq
                or "C"*bLen in bc4fq or "G"*bLen in bc4fq or "T"*bLen)):
            continue
        if(int(GetDescTagValue(fqRec.description, "FM")) < 2):
            continue
        if(workingBarcode == ""):
            workingBarcode = bc4fq
            workingSet = []
            workingSet.append(fqRec)
            continue
        elif(workingBarcode == bc4fq):
            workingSet.append(fqRec)
            continue
        elif(workingBarcode != bc4fq):
            mergedRecord, success = compareFastqRecords(
                workingSet, stringency=float(stringency))
            if(success):
                SeqIO.write(mergedRecord, outputHandle, "fastq")
            workingSet = []
            workingSet.append(fqRec)
            workingBarcode = bc4fq
            continue
        else:
            raise RuntimeError(
                "No idea what's going on. This code should be unreachable")
    inFq.close()
    outputHandle.close()
    return outFq


def TrimHoming(
        fq,
        homing,
        trimfq="default",
        bar_len=12,
        tags_file="default",
        trim_err="default",
        start_trim=1):
    pl("TrimHoming: \"{}\" from {}".format(homing, fq))
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
    if(trim_err == "default"):
        temp = '.'.join(fq.split('.')[0:-1] + ['err', 'fastq'])
        trim_err = temp.split('/')[-1]
    if(trimfq == "default"):
        temp = '.'.join(fq.split('.')[0:-1] + ['trim', 'fastq'])
        trimfq = temp.split('/')[-1]
    if(tags_file == "default"):
        temp = '.'.join(fq.split('.')[0:-1] + ['tags', 'fastq'])
        tags_file = temp.split('/')[-1]
    tagsOpen = open(tags_file, "w", 0)
    trimOpen = open(trimfq, "w", 0)
    InFastq = SeqIO.parse(fq, "fastq")
    HomingLen = len(homing)
    TotalTrim = HomingLen + bar_len + start_trim
    pl("Homing Length is {}".format(HomingLen))
    for rec in InFastq:
        pre_tag = SeqRecord(
            Seq(str(rec.seq)[0:bar_len], "fastq"),
            id=rec.id, description=rec.description)
        pre_tag.letter_annotations['phred_quality'] = rec.letter_annotations[
            'phred_quality'][0:bar_len]
        """
        if homing not in pre_tag.seq:
            pl("Homing sequence not in tag. Writing to error file.")
            SeqIO.write(rec,errOpen,"fastq")
            continue
        """
        SeqIO.write(pre_tag, tagsOpen, "fastq")
        descString = rec.description + " #G~BS=" + str(rec.seq[0:bar_len])
        post_tag = SeqRecord(Seq(str(rec.seq)[TotalTrim:],
                                 "fastq"),
                             id=rec.id,
                             description=descString)
        post_tag.letter_annotations['phred_quality'] = rec.letter_annotations[
            'phred_quality'][TotalTrim:]
        SeqIO.write(post_tag, trimOpen, "fastq")
    tagsOpen.close()
    trimOpen.close()
    return(tags_file, trimfq)
