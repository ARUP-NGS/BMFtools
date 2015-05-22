# cython: boundscheck=False, c_string_type=str, c_string_encoding=ascii
# cython: cdivision=True, cdivision_warnings=True, profile=True
from __future__ import division

"""
Contains various utilities for working with barcoded fastq files.
"""

import logging
import os
import shlex
import subprocess
from copy import copy as ccopy
import gzip
import sys
import collections
import time
import cStringIO
import operator
from operator import (add as oadd, le as ole, ge as oge, div as odiv,
                      mul as omul, add as oadd, attrgetter as oag)
from subprocess import check_call

import cython
import numpy as np
from numpy import (sum as nsum, amax as npamax, argmax as npargmax,
                   multiply as nmultiply, subtract as nsubtract, any as npany,
                   vstack as npvstack, greater as ngreater, less as nless,
                   array as nparray)
import pysam
from cytoolz import map as cmap, memoize
from pysam import fromQualityString

from utilBMF.HTSUtils import (PipedShellCall, GetSliceFastqProxy, ph2chr,
                              ph2chrDict, chr2ph, printlog as pl, FacePalm,
                              pFastqProxy, TrimExt)
from utilBMF import HTSUtils
from utilBMF.ErrorHandling import ThisIsMadness

cimport cython
cimport numpy as np
cimport pysam.cfaidx
cimport pysam.calignmentfile
cimport utilBMF.HTSUtils
ctypedef utilBMF.HTSUtils.pFastqProxy pFq


letterNumDict = {'A': 0, 'C': 0, 'G': 0, 'T': 0,
                 0: 'A', 1: 'C', 2: 'G', 3: 'T'}


@memoize
def chr2phFunc(x):
    return chr2ph[x]


@cython.locals(checks=cython.long, highMem=cython.bint,
               parallel=cython.bint, inFq1=cython.str,
               inFq2=cython.str, sortMem=cython.str)
def BarcodeSortBoth(inFq1, inFq2, sortMem="6G", parallel=False):
    if(parallel is False):
        pl("Parallel barcode sorting is set to false. Performing serially.")
        return BarcodeSort(inFq1), BarcodeSort(inFq2)
    outFq1 = '.'.join(inFq1.split('.')[0:-1] + ["BS", "fastq"])
    outFq2 = '.'.join(inFq2.split('.')[0:-1] + ["BS", "fastq"])
    pl("Sorting {} and {} by barcode sequence.".format(inFq1, inFq2))
    highMemStr = "-S " + sortMem
    BSstring1 = getBarcodeSortStr(inFq1, outFastq=outFq1,
                                  highMem=("6G" in sortMem))
    BSstring2 = getBarcodeSortStr(inFq2, outFastq=outFq2,
                                  highMem=("6G" in sortMem))
    pl("Background calling barcode sorting "
       "for read 1. Command: {}".format(BSstring1))
    BSCall1 = subprocess.Popen(BSstring1, stderr=None, shell=True,
                               stdout=None, stdin=None, close_fds=True)
    PipedShellCall(BSstring2)
    BSCall1.poll()
    checks = 0
    if(BSCall1.returncode == 0):
        return outFq1, outFq2
    elif(BSCall1.returncode is not None):
        raise subprocess.CalledProcessError("Barcode sort failed for read 1."
                                            " Call: {}".format(BSstring1))
    else:
        while checks < 300:
            checks = 0
            time.sleep(1)
            BSCall1.poll()
            if(BSCall1.returncode == 0):
                return outFq1, outFq2
            else:
                checks += 1
        raise subprocess.CalledProcessError("Barcode first sort didn't work, "
                                            "it seems - took more than 5 minu"
                                            "tes longer than the second barco"
                                            "de sort.")


@cython.locals(highMem=cython.bint)
def BarcodeSort(cython.str inFastq, cython.str outFastq="default",
                cython.bint highMem=True):
    cdef cython.str BSstring
    pl("Sorting {} by barcode sequence.".format(inFastq))
    BSstring = getBarcodeSortStr(inFastq, outFastq=outFastq, highMem=highMem)
    PipedShellCall(BSstring)
    pl("Barcode Sort shell call: {}".format(BSstring))
    if(outFastq == "default"):  # Added for compatibility with getBSstr
        outFastq = '.'.join(inFastq.split('.')[0:-1] + ["BS", "fastq"])
    return outFastq


@cython.returns(cython.str)
def getBarcodeSortStr(inFastq, outFastq="default", highMem=True):
    if(highMem):
        memStr = " -S 6G "
    if(outFastq == "default"):
        outFastq = '.'.join(inFastq.split('.')[0:-1] + ["BS", "fastq"])
    if(inFastq.endswith(".gz")):
        return ("zcat %s | paste - - - - | sort -t'|' -k3,3 -k1,1" % inFastq +
                " %s | tr '\t' '\n' > %s" % (memStr, outFastq))
    else:
        return ("cat %s | paste - - - - | sort -t'|' -k3,3 -k1,1" % inFastq +
                " %s | tr '\t' '\n' > %s" % (memStr, outFastq))


@cython.locals(stringency=cython.float, hybrid=cython.bint,
               famLimit=cython.long, keepFails=cython.bint,
               Success=cython.bint, PASS=cython.bint, frac=cython.float,
               compressB64=cython.bint, lenR=cython.long,
               numEq=cython.long, maxScore=cython.long, ND=cython.long)
def compareFqRecsFqPrx(list R, stringency=0.9, hybrid=False,
                       famLimit=1000, keepFails=True,
                       makeFA=True, makePV=True, name=None):
    """
    Compares the fastq records to create a consensus sequence (if it
    passes a filter)
    """
    if name is None:
        name = R[0].name
    cdef np.ndarray[np.int64_t, ndim = 1] phredQuals
    cdef np.ndarray[np.int64_t, ndim = 1] FA
    cdef np.ndarray[char, ndim = 1, mode="c"] finalSeq
    cdef cython.str lenRStr
    cdef list seqs
    cdef cython.str seq
    cdef cython.str seqItem
    cdef cython.str PVString
    cdef cython.str QualString
    cdef cython.str TagString
    cdef cython.str consFqString
    compress85 = True
    lenR = len(R)
    lenRStr = str(lenR)
    if(lenR > famLimit):
        logging.debug(
            "Read family - {} with {} members was capped at {}. ".format(
                R[0], lenR, famLimit))
        R = R[:famLimit]
    seqs = list(cmap(oag("sequence"), R))
    maxScore = 0
    Success = False
    numEq = 0
    for seq in seqs:
        numEq = sum(seq == seqItem for seqItem in seqs)
        if(oge(numEq, maxScore)):
            maxScore = numEq
            finalSeq = np.array(list(seq))
    try:
        frac = odiv(omul(numEq, 1.0), lenR)
    except ZeroDivisionError:
        pl("Length of R: {}".format(lenR))
        pl("numEq: {}".format(numEq))
    if(oge(frac, stringency)):
        Success = True
    elif(ole(frac, 0.5)):
        Success = False
    elif(hybrid):
        return compareFqRecsFast(R, makePV=makePV, makeFA=makeFA, name=name)
    FA = nparray([sum([seq[i] == finalSeq[i] for seq in seqs]) for i in
                  xrange(len(finalSeq))], dtype=np.int64)
    phredQuals = nparray([chr2ph[i] for i in list(R[0].quality)],
                         dtype=np.int64)
    phredQuals[phredQuals < 3] = 0
    finalSeq[phredQuals < 3] = "N"  # Let's see if this works!
    phredQuals = nmultiply(lenR, phredQuals, dtype=np.int64)
    if(npany(ngreater(phredQuals, 93))):
        QualString = "".join(list(cmap(ph2chr, phredQuals)))
        PVString = "|PV=" + ",".join(phredQuals.astype(str).tolist())
    else:
        QualString = "".join([ph2chrDict[i] for i in phredQuals])
        PVString = oadd("|PV=",
                        ",".join(phredQuals.astype(str).tolist()))
    TagString = "".join(["|FM=", lenRStr, "|FA=",
                         ",".join(nparray(FA).astype(str)),
                         "|ND=", str(nsubtract(lenR * len(seqs[0]),
                                                  nsum(FA))),
                         PVString])
    try:
        consFqString = "\n".join(
            ["".join(["@", name, " ", R[0].comment, TagString]),
             finalSeq.tostring(),
             "+",
             QualString, ""])
    except TypeError:
        print("TagString: {}".format(TagString))
        print("finalSeq: {}".format(finalSeq))
        print("QualString: {}".format(QualString))
        print("Name {}".format(R[0].name))
        print("Comment {}".format(R[0].comment))
        print("".join(["@", R[0].name, " ", R[0].comment, TagString]))
        raise ThisIsMadness("I can't figure out what's going on.")
    if(Success is False):
        return consFqString.replace("Pass", "Fail"), name
    return consFqString, name


@cython.returns(tuple)
def compareFqRecsFast(R, makePV=True, makeFA=True, name=None):
    """
    Calculates the most likely nucleotide
    at each position and returns the joined record string.
    """
    cdef cython.long lenR, ND
    cdef cython.bint Success
    cdef np.ndarray[np.int64_t, ndim = 2] quals
    cdef np.ndarray[np.int64_t, ndim = 2] qualA
    cdef np.ndarray[np.int64_t, ndim = 2] qualC
    cdef np.ndarray[np.int64_t, ndim = 2] qualG
    cdef np.ndarray[np.int64_t, ndim = 2] qualT
    cdef np.ndarray[np.int64_t, ndim = 2] qualAllSum
    cdef np.ndarray[np.int64_t, ndim = 1] MaxPhredSum
    cdef np.ndarray[np.int64_t, ndim = 1] phredQuals
    cdef np.ndarray[np.int64_t, ndim = 1] FA
    cdef np.ndarray[char, ndim = 1, mode="c"] newSeq
    cdef np.ndarray[cython.str, ndim = 1] seqs
    lenR = len(R)
    Success = True
    seqs = nparray(map(oag("sequence"), R), dtype=np.str_)
    stackArrays = tuple([np.char.array(s, itemsize=1) for s in seqs])
    seqArray = npvstack(stackArrays)
    if(name is None):
        name = R[0]

    # print(repr(seqArray))
    quals = nparray(map(fromQualityString, map(oag("quality"), R)))
    """
    quals = nparray(
        [list(cmap(chr2phFunc, list(record.quality))) for record in R],
        dtype=np.int64)
    """
    # Qualities of 2 are placeholders and mean nothing in Illumina sequencing.
    # Let's turn them into what they should be: nothing.
    quals[quals < 3] = 0
    qualA = ccopy(quals)
    qualC = ccopy(quals)
    qualG = ccopy(quals)
    qualT = ccopy(quals)
    qualA[seqArray != "A"] = 0
    qualA = nsum(qualA, 0, dtype=np.int64)
    qualC[seqArray != "C"] = 0
    qualC = nsum(qualC, 0, dtype=np.int64)
    qualG[seqArray != "G"] = 0
    qualG = nsum(qualG, 0, dtype=np.int64)
    qualT[seqArray != "T"] = 0
    qualT = nsum(qualT, 0, dtype=np.int64)
    qualAllSum = npvstack(
        [qualA, qualC, qualG, qualT])
    newSeq = np.array([letterNumDict[i] for i in npargmax(qualAllSum, 0)])
    MaxPhredSum = npamax(qualAllSum, 0)  # Avoid calculating twice.
    phredQuals = nsubtract(nmultiply(2, MaxPhredSum, dtype=np.int64),
                           nsum(qualAllSum, 0, dtype=np.int64),
                           dtype=np.int64)
    newSeq[phredQuals == 0] = "N"
    FA = nparray([sum([seq[i] == newSeq[i] for seq in seqs]) for i in
                  xrange(len(newSeq))], dtype=np.int64)
    ND = lenR * len(seqs[0]) - nsum(FA)
    if(npany(nless(phredQuals, 0))):
        pl("repr of phredQuals %s" % repr(phredQuals), level=logging.DEBUG)
        phredQuals = abs(phredQuals)
    if(npany(ngreater(phredQuals, 93))):
        PVString = "|PV=" + ",".join(phredQuals.astype(str))
        phredQuals[phredQuals > 93] = 93
        phredQualsStr = "".join([ph2chrDict[i] for i in phredQuals])
    else:
        phredQualsStr = "".join([ph2chrDict[i] for i in phredQuals])
        PVString = "|PV=" + ",".join(phredQuals.astype(str))
    TagString = "".join(["|FM=", str(lenR), "|ND=",
                         str(ND), "|FA=", FA.astype(str), PVString])
    consolidatedFqStr = "\n".join([
        "".join(["@", name, " ", R[0].comment, TagString]),
        newSeq.tostring(),
        "+",
        phredQualsStr, ""])
    if(not Success):
        return consolidatedFqStr.replace("Pass", "Fail"), name
    return consolidatedFqStr, name


@cython.returns(cython.str)
def CutadaptPaired(cython.str fq1, cython.str fq2,
                   p3Seq="default", p5Seq="default",
                   cython.long overlapLen=6, cython.bint makeCall=True):
    """
    Returns a string which can be called for running cutadapt v.1.7.1
    for paired-end reads in a single call.
    """
    cdef cython.str outfq1
    cdef cython.str outfq2
    cdef cython.str commandStr
    outfq1 = ".".join(fq1.split('.')[0:-1] + ["cutadapt", "fastq"])
    outfq2 = ".".join(fq2.split('.')[0:-1] + ["cutadapt", "fastq"])
    if(p3Seq == "default"):
        HTSUtils.FacePalm("3-prime primer sequence required for cutadapt!")
    if(p5Seq == "default"):
        pl("No 5' sequence provided for cutadapt. Only trimming 3'.")
        commandStr = ("cutadapt --mask-adapter --match-read-wildcards"
                      "-a {} -o {} -p {} -O {} {} {}".format(p3Seq,
                                                             outfq1,
                                                             outfq2,
                                                             overlapLen,
                                                             fq1, fq2))
    else:
        commandStr = ("cutadapt --mask-adapter --match-read-wildcards -a "
                      "{} -g {} -o {} -p".format(p3Seq, p5Seq, outfq1) +
                      " {} -O {} {} {}".format(outfq2, overlapLen, fq1, fq2))
    pl("Cutadapt command string: {}".format(commandStr))
    if(makeCall):
        subprocess.check_call(shlex.split(commandStr))
        return outfq1, outfq2
    return commandStr


@cython.locals(overlapLen=cython.long)
@cython.returns(cython.str)
def CutadaptString(fq, p3Seq="default", p5Seq="default", overlapLen=6):
    """
    Returns a string which can be called for running cutadapt v.1.7.1.
    """
    outfq = ".".join(fq.split('.')[0:-1] + ["cutadapt", "fastq"])
    if(p3Seq == "default"):
        HTSUtils.FacePalm("3-prime primer sequence required for cutadapt!")
    if(p5Seq == "default"):
        pl("No 5' sequence provided for cutadapt. Only trimming 3'.")
        commandStr = "cutadapt --mask-adapter {} -a {} -o {} -O {} {}".format(
            "--match-read-wildcards", p3Seq, outfq, overlapLen, fq)
    else:
        commandStr = ("cutadapt --mask-adapter --match-read-wildcards -a "
                      "{} -g {} -o {} -O {} {}".format(p3Seq, p5Seq, outfq,
                                                       overlapLen, fq))
    pl("Cutadapt command string: {}".format(commandStr))
    return commandStr, outfq


@cython.locals(overlapLen=cython.long)
def CallCutadapt(fq, p3Seq="default", p5Seq="default", overlapLen=6):
    """
    Calls cutadapt to remove adapter sequence at either end of the reads.
    Written for v1.7.1 and single-end calls.
    """
    commandStr, outfq = CutadaptString(fq, p3Seq=p3Seq, p5Seq=p5Seq,
                                       overlapLen=overlapLen)
    subprocess.check_call(commandStr)
    return outfq


@cython.locals(overlapLap=cython.long, numChecks=cython.long)
def CallCutadaptBoth(fq1, fq2, p3Seq="default", p5Seq="default", overlapLen=6):
    fq1Str, outfq1 = CutadaptString(fq1, p3Seq=p3Seq, p5Seq=p5Seq,
                                    overlapLen=overlapLen)
    fq2Str, outfq2 = CutadaptString(fq2, p3Seq=p3Seq, p5Seq=p5Seq,
                                    overlapLen=overlapLen)
    pl("About to open a Popen instance for read 2. Command: {}".format(fq2Str))
    fq2Popen = subprocess.Popen(fq2Str, shell=True, stdin=None, stderr=None,
                                stdout=None, close_fds=True)
    pl("Cutadapt running in the background for read 2. Now calling for read "
       "1: {}".format(fq1Str))
    subprocess.check_call(fq1Str, shell=True)
    numChecks = 0
    while True:
        if(numChecks >= 300):
            raise subprocess.CalledProcessError(
                "Cutadapt took more than 5 minutes longer for read 2 than rea"
                "d 1. Did something go wrong?")
        if fq2Popen.poll() == 0:
            return outfq1, outfq2
        elif(fq2Popen.poll() is None):
            pl("Checking if cutadapt for read 2 is finished. Seconds elapsed ("
               "approximate): {}".format(numChecks))
            numChecks += 1
            time.sleep(1)
            continue
        else:
            raise subprocess.CalledProcessError("Cutadapt failed for read 2!")


@cython.locals(useGzip=cython.bint, bLen=cython.long)
def FastqPairedShading(fq1, fq2, indexfq="default",
                       useGzip=False, readPairsPerWrite=10,
                       cython.long head=2):
    """
    Tags fastqs with barcodes from an index fastq.
    """
    #  C declarations
    cdef pysam.cfaidx.FastqProxy read1
    cdef pysam.cfaidx.FastqProxy read2
    cdef pysam.cfaidx.FastqProxy indexRead
    cdef cython.str outfq1
    cdef cython.str outfq2
    pl("Now beginning fastq marking: Pass/Fail and Barcode")
    if(indexfq == "default"):
        raise ValueError("For an i5/i7 index ")
    outfq1 = TrimExt(fq1).replace(".fastq",
                                  "").split("/")[-1] + ".shaded.fastq"
    outfq2 = TrimExt(fq2).replace(".fastq",
                                  "").split("/")[-1] + ".shaded.fastq"
    if(useGzip):
        outfq1 += ".gz"
        outfq2 += ".gz"
    pl("Output fastqs: {}, {}.".format(outfq1, outfq2))
    inFq1 = pysam.FastqFile(fq1)
    inFq2 = pysam.FastqFile(fq2)
    ifn1 = inFq1.next
    ifn2 = inFq2.next
    if useGzip is False:
        outFqHandle1 = open(outfq1, "w")
        outFqHandle2 = open(outfq2, "w")
        f1 = cStringIO.StringIO()
        f2 = cStringIO.StringIO()
    else:
        outFqHandle1 = open(outfq1, "wb")
        outFqHandle2 = open(outfq2, "wb")
        cString1 = cStringIO.StringIO()
        cString2 = cStringIO.StringIO()
        f1 = gzip.GzipFile(fileobj=cString1, mode="w")
        f2 = gzip.GzipFile(fileobj=cString2, mode="w")
    inIndex = pysam.FastqFile(indexfq)
    ifin = inIndex.next
    outFqSet1 = []
    outFqSet2 = []
    numWritten = 0
    ofh1w = outFqHandle1.write
    ofh2w = outFqHandle2.write
    while True:
        if(numWritten >= readPairsPerWrite):
            if(not useGzip):
                ofh1w(f1.getvalue())
                ofh2w(f2.getvalue())
                f1 = cStringIO.StringIO()
                f2 = cStringIO.StringIO()
            else:
                f1.flush()
                f2.flush()
                ofh1w(cString1.getvalue())
                ofh2w(cString2.getvalue())
                cString1 = cStringIO.StringIO()
                cString2 = cStringIO.StringIO()
                f1 = gzip.GzipFile(fileobj=cString1, mode="w")
                f2 = gzip.GzipFile(fileobj=cString2, mode="w")
            numWritten = 0
        try:
            read1 = ifn1()
        except StopIteration:
            break
        read2 = ifn2()
        indexRead = ifin()
        tempBar = indexRead.sequence
        bLen = len(indexRead.sequence) * 5 // 6
        # bLen - 10 of 12 in a row, or 5/6. See Loeb, et al.
        # This is for removing low complexity reads
        # print("bLen is {}".format(bLen))
        if(("N" in tempBar or "A" * bLen in tempBar or
                "C" * bLen in tempBar or "G" * bLen in tempBar or
                "T" * bLen in tempBar)):
            '''
            pl("Failing barcode for read {} is {} ".format(indexRead,
                                                           tempBar),
               level=logging.DEBUG)
            '''
            tempBar = "%s%s%s" % (read1.sequence[:head], tempBar, read2.sequence[:head])
            f1.write("@%s %s|FP=IndexFail|BS=" % (read1.name, read1.comment) +
                     "%s\n%s\n+\n%s\n" % (tempBar, read1.sequence, read1.quality))
            f2.write("@%s %s|FP=IndexFail|BS=" % (read1.name, read2.comment) +
                     "%s\n%s\n+\n%s\n" % (tempBar, read2.sequence, read2.quality))
        else:
            tempBar = "%s%s%s" % (read1.sequence[:head], tempBar, read2.sequence[:head])
            f1.write("@%s %s|FP=IndexPass|BS=" % (read1.name, read1.comment) +
                     "%s\n%s\n+\n%s\n" % (tempBar, read1.sequence, read1.quality))
            f2.write("@%s %s|FP=IndexPass|BS=" % (read1.name, read2.comment) +
                     "%s\n%s\n+\n%s\n" % (tempBar, read2.sequence, read2.quality))
        numWritten += 1
    if(useGzip is False):
        ofh1w(f1.getvalue())
        ofh2w(f2.getvalue())
    else:
        f1.close()
        f2.close()
        ofh1w(cString1.getvalue())
        ofh2w(cString2.getvalue())
    outFqHandle1.close()
    outFqHandle2.close()
    return outfq1, outfq2


def FastqSingleShading(fq,
                       indexfq="default",
                       outfq="default",
                       cython.bint gzip=False,
                       cython.long head=0):
    cdef pysam.cfaidx.FastqProxy read1
    cdef pFq pRead1, pIndexRead
    cdef pysam.cfaidx.FastqFile inFq1, inIndex
    pl("Now beginning fastq marking: Pass/Fail and Barcode")
    if(indexfq == "default"):
        raise ValueError("For an i5/i7 index ")
    if(outfq == "default"):
        outfq = '.'.join(fq.split('.')[0:-1]) + '.shaded.fastq'
    inFq1 = pysam.FastqFile(fq)
    outFqHandle1 = open(outfq, "w")
    inIndex = pysam.FastqFile(indexfq)
    for read1 in inFq1:
        pIndexRead = pFastqProxy(inIndex.next())
        pRead1 = pFastqProxy(read1)
        if("N" in pRead1.pIndexRead.sequence):
            read1.comment += "|FP=IndexFail|BS=%" % (read1.sequence[:head] +
                                                     pIndexRead.sequence)
            outFqHandle1.write(str(read1))
        else:
            read1.comment += "|FP=IndexPass|BS=%s" % (read1.sequence[:head] +
                                                      pIndexRead.sequence)
            outFqHandle1.write(str(read1))
    outFqHandle1.close()
    if(gzip):
        check_call(['gzip', fq], shell=False)
    return


def HomingSeqLoc(fq, homing, bcLen=12):
    pl("Now beginning HomingSeqLoc.")
    cdef pysam.cfaidx.FastqProxy read
    cdef utilBMF.HTSUtils.pFastqProxy pRead
    InFastq = pysam.FastqFile(fq)
    Tpref = '.'.join(fq.split('.')[0:-1])
    Prefix = Tpref.split('/')[-1]
    StdFilename = Prefix + '.{}.fastq'.format("homing" + str(bcLen))
    ElseFilename = Prefix + '.else.fastq'
    ElseLoc = Prefix + '.else.supp'
    StdFastq = open(StdFilename, 'w', 0)  # Homing at expected Location
    ElseFastq = open(ElseFilename, 'w', 0)
    ElseLocations = open(ElseLoc, 'w', 0)
    for read in InFastq:
        pRead = pFastqProxy(read)
        seq = pRead.sequence
        if(seq.find(homing) == -1):
            pRead.comment += "|FP=HomingFail"
            ElseFastq.write(str(pRead))
        elif(seq[bcLen:bcLen + len(homing)] == homing):
            pRead.comment += "|FP=HomingPass"
            StdFastq.write(str(pRead))
        else:
            pRead.comment = "|FP=HomingFail"
            ElseFastq.write(str(pRead))
            ElseLocations.write(repr(seq.find(homing)) + "\t" +
                                pRead.name + "\n")
    StdFastq.close()
    ElseFastq.close()
    ElseLocations.close()
    return StdFilename, ElseFilename


def fastq_sort(in_fastq, out_fastq):
    pl("Now beginning fastq_sort.")
    outfile = open(out_fastq, 'w')
    command_str = ('cat ' + in_fastq + ' | paste - - - - | '
                   'sort -k1,1 -t " " | tr "\t" "\n"')
    subprocess.check_call(command_str, stdout=open(outfile, "w"), shell=True)
    outfile.close()
    return(command_str)


def FastqRegex(fq, string, matchFile="default", missFile="default"):
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


@cython.returns(cython.str)
def GetDescTagValue(readDesc, tag="default"):
    """
    Gets the value associated with a given tag.
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
        raise KeyError("Invalid tag: %s" % tag)


def GetDescriptionTagDict(readDesc):
    """Returns a set of key/value pairs in a dictionary for """
    tagSetEntries = [i.strip().split("=") for i in readDesc.split("|")][1:]
    tagDict = {}
    try:
        for pair in tagSetEntries:
            tagDict[pair[0]] = pair[1].split(' ')[0]
    except IndexError:
        pl("A value is stored with the| tag which doesn't contain an =.")
        pl("tagSetEntries: {}".format(tagSetEntries))
        raise IndexError("Check that fastq description meets specifications.")
    # pl("Repr of tagDict is {}".format(tagDict))
    except TypeError:
        pl("tagSetEntries: {}".format(tagSetEntries))
        raise ThisIsMadness("YOU HAVE NO CHANCE TO SURVIVE MAKE YOUR TIME")
    return tagDict


def pairedFastqConsolidate(fq1, fq2, cython.float stringency=0.9,
                           cython.long readPairsPerWrite=100,
                           cython.bint UsecProfile=False,
                           cython.bint onlyNumpy=True,
                           cython.bint skipSingles=False,
                           cython.bint skipFails=False):
    if(UsecProfile):
        import cProfile
        import pstats
        pr = cProfile.Profile()
        pr.enable()
    cdef cython.str outFqPair1, outFqPair2, workingBarcode, bc4fq1
    cdef pysam.cfaidx.FastqFile inFq1, inFq2
    cdef pFq fqRec, fqRec2
    cdef list workingSet1, workingSet2, StringList1, StringList2
    cdef cython.long numProc
    outFqPair1 = TrimExt(fq1) + ".cons.fastq"
    outFqPair2 = TrimExt(fq2) + '.cons.fastq'
    pl("Now running pairedFastqConsolidate on {} and {}.".format(fq1,fq2))
    pl("Command required to duplicate this action:"
       " pairedFastqConsolidate('{}', '{}', ".format(fq1, fq2) +
       "stringency={}, readPairsPerWrite={})".format(stringency,
                                                     readPairsPerWrite))
    inFq1 = pysam.FastqFile(fq1)
    inFq2 = pysam.FastqFile(fq2)
    outputHandle1 = open(outFqPair1, 'w')
    outputHandle2 = open(outFqPair2, 'w')
    # cString1 = cStringIO.StringIO()
    # cString2 = cStringIO.StringIO()
    StringList1 = []
    StringList2 = []
    workingBarcode = ""
    workingSet1 = []
    workingSet2 = []
    numProc = 0
    while True:
        if(numProc % readPairsPerWrite == 0):
            # outputHandle1.write(cString1.getvalue())
            # outputHandle2.write(cString2.getvalue())
            outputHandle1.write("".join(StringList1))
            outputHandle2.write("".join(StringList2))
            StringList1 = []
            StringList2 = []
            sl1a = StringList1.append
            sl2a = StringList2.append
            # cString1 = cStringIO.StringIO()
            # cString2 = cStringIO.StringIO()
        try:
            fqRec = pFastqProxy(inFq1.next())
        except StopIteration:
            break
        bc4fq1 = GetDescTagValue(fqRec.comment, "BS")
        fqRec2 = pFastqProxy(inFq2.next())
        # Originally removing reads with family size <2, since one pair could
        # have more than the other, it's important that I keep these reads in
        # and filter them from the BAM file
        if(workingBarcode == ""):
            try:
                workingBarcode = bc4fq1
                workingSet1 = [fqRec]
                workingSet2 = [fqRec2]
                continue
            except TypeError:
                print("workingBarcode = " + workingBarcode)
                print("workingSet1 = " + repr(workingSet1))
                print("workingSet2 = " + repr(workingSet2))
                sys.exit()
        elif(workingBarcode == bc4fq1):
            workingSet1.append(fqRec)
            workingSet2.append(fqRec2)
            continue
        elif(workingBarcode != bc4fq1):
            if(skipSingles and len(workingSet1) == 1):
                workingBarcode = ""
                workingSet1 = []
                workingSet2 = []
                continue
            # cString1.write(compareFqRecsFqPrx(workingSet1) + "\n")
            # cString2.write(compareFqRecsFqPrx(workingSet2) + "\n")
            # String1 += compareFqRecsFqPrx(workingSet1) + "\n"
            # String2 += compareFqRecsFqPrx(workingSet2) + "\n"
            tStr1, name = compareFqRecsFqPrx(workingSet1)
            tStr2, name = compareFqRecsFqPrx(workingSet2, name=name)
            if(skipFails and ("Fail" in tStr1 or "Fail" in tStr2)):
                continue
            sl1a(tStr1)
            sl2a(tStr2)
            workingSet1 = [fqRec]
            workingSet2 = [fqRec2]
            workingBarcode = bc4fq1
            numProc += 1
            continue
    if(UsecProfile):
        s = cStringIO.StringIO()
        pr.disable()
        sortby = "cumulative"
        ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
        ps.print_stats()
        open("cProfile.stats.txt", "w").write(s.getvalue())
    # outputHandle1.write(cString1.getvalue())
    # outputHandle2.write(cString2.getvalue())
    outputHandle1.write("".join(StringList1))
    outputHandle2.write("".join(StringList2))
    outputHandle1.flush()
    outputHandle2.flush()
    inFq1.close()
    inFq2.close()
    # cString1.close()
    # cString2.close()
    outputHandle1.close()
    outputHandle2.close()
    return outFqPair1, outFqPair2


def TrimHomingSingle(
        fq,
        cython.str homing=None,
        cython.str trimfq=None,
        cython.long bcLen=12,
        cython.str trim_err=None,
        cython.long start_trim=1):
    cdef pysam.cfaidx.FastqProxy read
    cdef cython.long HomingLen, TotalTrim
    pl("TrimHoming: \"{}\" from {}.".format(homing, fq))
    if(trim_err is None):
        trim_err = TrimExt(fq) + '.err.fastq'
    if(trimfq is None):
        trimfq = TrimExt(fq) + ".trim.fastq"
    trimHandle = open(trimfq, "w")
    errHandle = open(trim_err, "w")
    InFastq = pysam.FastqFile(fq)
    HomingLen = len(homing)
    TotalTrim = HomingLen + bcLen + start_trim
    tw = trimHandle.write
    ew = errHandle.write
    for read in InFastq:
        homingLoc = read.sequence[bcLen:bcLen + HomingLen]
        if homing not in homingLoc:
            pl("Homing sequence not in tag. Writing to error file.",
               level=logging.DEBUG)
            ew(str(pFastqProxy(read)))
            continue
        tw(GetSliceFastqProxy(read, firstBase=TotalTrim,
                              addString="|BS=" + read.sequence[0:bcLen]))
    trimHandle.close()
    errHandle.close()
    return trimfq


def TrimHomingPaired(inFq1, inFq2, cython.long bcLen=12,
                     cython.str homing=None, cython.str trimfq1=None,
                     cython.str trimfq2=None, cython.long start_trim=1):
    cdef pysam.cfaidx.FastqProxy read1, read2
    cdef cython.long HomingLen, TotalTrim
    pl("Getting inline barcodes for files %s, %s with homing %s." % (inFq1,
                                                                     inFq2,
                                                                     homing))
    trim_err = TrimExt(inFq1) + '.err.fastq'
    if(trimfq1 is None):
        trimfq1 = TrimExt(inFq1) + ".trim.fastq"
    if(trimfq2 is None):
        trimfq2 = TrimExt(inFq2) + ".trim.fastq"
    trimHandle1 = open(trimfq1, "w")
    trimHandle2 = open(trimfq2, "w")
    errHandle = open(trim_err, "w")
    InFastq1 = pysam.FastqFile(inFq1)
    InFastq2 = pysam.FastqFile(inFq2)
    fqNext = InFastq1.next
    HomingLen = len(homing)
    TotalTrim = HomingLen + bcLen + start_trim
    pl("Homing length: %s. TotalTrim: %s" % (HomingLen, TotalTrim))
    tw1 = trimHandle1.write
    tw2 = trimHandle2.write
    ew = errHandle.write
    for read1 in InFastq1:
        read2 = fqNext()
        if homing not in read1.sequence[bcLen:bcLen + HomingLen]:
            ew(str(pFastqProxy(read1)))
            ew(str(pFastqProxy(read2)))
            continue
        if homing not in read2.sequence[bcLen:bcLen + HomingLen]:
            ew(str(pFastqProxy(read1)))
            ew(str(pFastqProxy(read2)))
            continue
        barcode = read1.sequence[0:bcLen] + read2.sequence[0:bcLen]
        tw1(GetSliceFastqProxy(read1, firstBase=TotalTrim,
                               addString="|BS=%s" % barcode))
        tw2(GetSliceFastqProxy(read2, firstBase=TotalTrim,
                               addString="|BS=%s" % barcode))
    trimHandle1.close()
    trimHandle2.close()
    errHandle.close()
    return trimfq1, trimfq2


@cython.locals(asDict=cython.bint)
@cython.returns(cython.str)
def CalcFamUtils(inFq, asDict=False):
    """
    Uses bioawk to get summary statistics on a fastq file quickly.
    """
    commandStr = ("bioawk -c fastx {{n=split($comment,array," ");k=array[4];n="
                  "split(k,array,\"=\"); len += 1; sum += array[2]; if(array[2"
                  "] == 1) {{singletons += 1}};if(array[2] >= 2) {{realFams +="
                  " 1; readFmSum += array[2]}}}};END {{print \"MFS=\"sum /len"
                  "\";NumRealFams=\"realFams\";ReadFamSum=\"readFmSum/realFams"
                  "\";NumSingletons=\"singletons}}' %s" % inFq)
    strOut = subprocess.check_output(shlex.split(commandStr)).strip()
    if(asDict):
        return dict([i.split("=") for i in strOut.split(";")])
    return strOut
