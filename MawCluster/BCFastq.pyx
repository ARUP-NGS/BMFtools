# cython: c_string_type=str, c_string_encoding=ascii
# cython: cdivision=True, cdivision_warnings=True, profile=True
from __future__ import division

"""
Contains various utilities for working with barcoded fastq files.
"""

import logging
import os
import shlex
import subprocess
from subprocess import check_output
from copy import copy as ccopy
import gzip
import sys
import collections
import time
import cStringIO
import operator
import uuid
from operator import (add as oadd, le as ole, ge as oge, div as odiv,
                      mul as omul, add as oadd, attrgetter as oag,
                      methodcaller as mc)
from subprocess import check_call

import cython
import numpy as np
from numpy import (sum as nsum, amax as npamax, argmax as npargmax,
                   multiply as nmul, subtract as nsub, any as npany,
                   vstack as npvstack, greater as ngreater, less as nless,
                   array as nparray, divide as ndiv, char as npchar)
import pysam
from cytoolz import memoize
from functools import partial
from itertools import groupby

from utilBMF.HTSUtils import (SliceFastqProxy, ph2chr,
                              printlog as pl,
                              pFastqProxy, TrimExt, pFastqFile, getBS,
                              hamming_cousins)
from utilBMF import HTSUtils
from utilBMF.ErrorHandling import ThisIsMadness as Tim, FunctionCallException
try:
    import re2 as re
except ImportError:
    import re
oagseq = oag("sequence")
oagqual = oag("quality")
npchararray = npchar.array
partialnpchar = partial(npchararray, itemsize=1)


letterNumDict = {'A': 0, 'C': 0, 'G': 0, 'T': 0,
                 0: 'A', 1: 'C', 2: 'G', 3: 'T'}


@memoize
def chr2phFunc(x):
    return chr2ph[x]


def SortAndMarkFastqsCommand(Fq1, Fq2, IndexFq):
    return ("pr -mts <(cat %s | paste - - - -) <(cat %s | " % (Fq1, Fq2) +
            "paste - - - -) <(cat %s | paste - - - -) | awk" % IndexFq +
            "'BEGIN {{FS=OFS=\"\t\"}};{{BS=\"|BS=\"$10\"\"substr($2,"
            "0,2)\"\"substr($6,0,2); print $1\"\"BS, $2, $3, $4, $5\"\"BS, "
            "$6, $7, $8}}' | sort -t\"|\" -k2,2 -k1,1 | tr '\t' '\n'")


@cython.locals(checks=int,
               parallel=cython.bint, sortMem=cystr)
def BarcodeSortBoth(cystr inFq1, cystr inFq2,
                    cystr sortMem="6G", cython.bint parallel=False):
    cdef cystr outFq1, outFq2, highMemStr
    if(parallel is False):
        pl("Parallel barcode sorting is set to false. Performing serially.")
        return BarcodeSort(inFq1), BarcodeSort(inFq2)
    outFq1 = '.'.join(inFq1.split('.')[0:-1] + ["BS", "fastq"])
    outFq2 = '.'.join(inFq2.split('.')[0:-1] + ["BS", "fastq"])
    pl("Sorting {} and {} by barcode sequence.".format(inFq1, inFq2))
    highMemStr = "-S " + sortMem
    BSstring1 = getBarcodeSortStr(inFq1, outFastq=outFq1,
                                  mem=sortMem)
    BSstring2 = getBarcodeSortStr(inFq2, outFastq=outFq2,
                                  mem=sortMem)
    pl("Background calling barcode sorting "
       "for read 1. Command: {}".format(BSstring1))
    BSCall1 = subprocess.Popen(BSstring1, stderr=None, shell=True,
                               stdout=None, stdin=None, close_fds=True)
    check_call(BSstring2, shell=True)
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
def BarcodeSort(cystr inFastq, cystr outFastq="default",
                cystr mem="6G"):
    cdef cystr BSstring
    pl("Sorting {} by barcode sequence.".format(inFastq))
    BSstring = getBarcodeSortStr(inFastq, outFastq=outFastq, mem=mem)
    check_call(BSstring, shell=True)
    pl("Barcode Sort shell call: {}".format(BSstring))
    if(outFastq == "default"):  # Added for compatibility with getBSstr
        outFastq = '.'.join(inFastq.split('.')[0:-1] + ["BS", "fastq"])
    return outFastq


@cython.returns(cystr)
def getBarcodeSortStr(inFastq, outFastq="default", mem=""):
    if(mem != ""):
        mem = " -S " + mem
    if(outFastq == "default"):
        outFastq = '.'.join(inFastq.split('.')[0:-1] + ["BS", "fastq"])
    if(inFastq.endswith(".gz")):
        return ("zcat %s | paste - - - - | sort -t'|' -k3,3 -k1,1" % inFastq +
                " %s | tr '\t' '\n' > %s" % (mem, outFastq))
    else:
        return ("cat %s | paste - - - - | sort -t'|' -k3,3 -k1,1" % inFastq +
                " %s | tr '\t' '\n' > %s" % (mem, outFastq))


@cython.boundscheck(False)
cpdef cystr cFRP_helper(list R, cystr name=None):
    return compareFqRecsFqPrx(R, name=name)


@cython.boundscheck(False)
cdef cystr compareFqRecsFqPrx(list R, cystr name=None,
                              float stringency=0.9,
                              int famLimit=1000,
                              cython.bint keepFails=True,
                              object oagseq=oagseq,
                              dict chr2ph=chr2ph,
                              dict ph2chrDict=ph2chrDict,
                              object ph2chr=ph2chr,
                              object int2Str=int2Str):
    """
    TODO: Unit test for this function.
    Compares the fastq records to create a consensus sequence (if it
    passes a filter).
    If NLowQual is set, all bases with q < 3 are set to N.
    If hybrid is set, a failure to successfully demultiplex falls back to a
    base by base comparison.
    """
    cdef ndarray[np.int64_t, ndim=1] phredQuals, FA
    cdef ndarray[char, ndim=1, mode = "c"] finalSeq
    cdef list seqs
    cdef cystr seqItem, seq, lenRStr, qual
    cdef cystr PVString, QualString, TagString, consFqString
    cdef int lenR, numEq, maxScore
    cdef cython.bint Success
    if(name is None):
        name = R[0].name
    lenR = len(R)
    lenRStr = str(lenR)
    '''
    if lenR == 1:
        PVString = "|PV=%s" % ",".join([chr2phStr[i] for i in R[0].quality])
        TagString = "|FM=1|ND=0|FA=%s|PV=%s" % (",".join(["1"] * lenSeq),
                                                PVString)
        return "@%s %s%s\n%s\n+\n%s\n" % (R[0].name, R[0].comment,
                                          TagString, R[0].sequence,
                                          R[0].quality)
    '''
    if(lenR > famLimit):
        logging.debug(
            "Read family - {} with {} members was capped at {}. ".format(
                R[0], lenR, famLimit))
        R = R[:famLimit]
    seqs = map(oagseq, R)
    maxScore = 0
    Success = False
    numEq = 0
    for seq in seqs:
        numEq = sum([seq == seqItem for seqItem in seqs])
        if(numEq > maxScore):
            maxScore = numEq
            finalSeq = np.array(list(seq))
    try:
        frac = numEq * 1. / lenR
    except ZeroDivisionError:
        pl("Length of R: {}".format(lenR))
        pl("numEq: {}".format(numEq))
    if(frac > stringency):
        Success = True
    elif(frac < 0.5):
        Success = False
    FA = nparray([sum([seq[i] == finalSeq[i] for seq in seqs]) for i in
                 xrange(len(finalSeq))], dtype=np.int64)
    phredQuals = nparray([sum([chr2ph[qual[i]] for
                               qual in map(oagqual, R) if
                               seq[i] == finalSeq[i]]) for
                          i in xrange(len(finalSeq))])
    # phredQuals[phredQuals < 3] = 0  # Throw array q scores of 2.
    # finalSeq[phredQuals < 3] = "N"  # Set all bases with q < 3 to N
    try:
        QualString = "".join([ph2chrDict[i] for i in phredQuals])
    except KeyError:
        QualString = "".join(map(ph2chr, phredQuals))
    try:
        PVString = "|PV=%s" % (",".join([int2Str[i] for i in phredQuals]))
    except KeyError:
        PVString = "|PV=%s" % (",".join(phredQuals.astype(str)))
    try:
        FAString = "|FA=%s" % (",".join([int2Str[i] for i in FA]))
    except KeyError:
        FAString = "|FA=%s" % (",".join(FA.astype(str)))
    TagString = "|FM=%s%s|ND=%s%s" % (
        lenRStr,  FAString, lenR * len(finalSeq) - nsum(FA), PVString)
    """
    try:
    """
    consFqString = "@%s %s%s\n%s\n+\n%s\n" % (name, R[0].comment,
                                              TagString,
                                              finalSeq.tostring(), QualString)
    """
    Uncomment and add indent the consFqString assignment command
    for debugging.
    except TypeError:
        print("TagString: {}".format(TagString))
        print("finalSeq: {}".format(finalSeq))
        print("QualString: {}".format(QualString))
        print("Name {}".format(R[0].name))
        print("Comment {}".format(R[0].comment))
        print("".join(["@", R[0].name, " ", R[0].comment, TagString]))
        raise Tim("I can't figure out what's going on.")
    """
    if(Success is False):
        return consFqString.replace("Pass", "Fail")
    return consFqString


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef cystr cFRF_helper(list R, cystr name=None):
    return compareFqRecsFast(R, name=name)


@cython.boundscheck(False)
@cython.wraparound(False)
cdef cystr compareFqRecsFast(list R,
                             cystr name=None,
                             int famLimit=100,
                             dict chr2ph=chr2ph,
                             dict letterNumDict=letterNumDict,
                             dict ph2chrDict=ph2chrDict,
                             object ccopy=ccopy,
                             object npchararray=npchararray,
                             object oagqual=oagqual,
                             object oagseq=oagseq,
                             object partialnpchar=partialnpchar,
                             object ph2chr=ph2chr,
                             object chr2phStr=chr2phStr,
                             object int2Str=int2Str):
    """
    TODO: Unit test for this function.
    Also, consider making a cpdef version!
    Calculates the most likely nucleotide
    at each position and returns the joined record string.
    """
    cdef int lenR, ND, lenSeq
    cdef cython.bint Success
    cdef cystr seq, qual, seqItem, qualChar, PVString, TagString
    cdef cystr consolidatedFqStr
    cdef ndarray[np.int64_t, ndim=2] quals, qualA, qualC, qualG
    cdef ndarray[np.int64_t, ndim=2] qualT, qualAllSum
    cdef ndarray[np.int64_t, ndim=1] qualAFlat, qualCFlat, qualGFlat, FA
    cdef ndarray[np.int64_t, ndim=1] MaxPhredSum, phredQuals, qualTFlat
    cdef ndarray[char, ndim=1, mode = "c"] newSeq
    if(name is None):
        name = R[0].name
    lenR = len(R)
    lenSeq = len(R[0].sequence)
    if lenR == 1:
        TagString = "|FM=1|ND=0|FA=%s|PV=%s" % (",".join(["1"] * lenSeq),
                                                ",".join([chr2phStr[i] for
                                                          i in R[0].quality]))
        return "@%s %s%s\n%s\n+\n%s\n" % (name, R[0].comment,
                                          TagString, R[0].sequence,
                                          R[0].quality)
    elif lenR > famLimit:
        return compareFqRecsFqPrx(R, name=name)  # Lazy large fams
    Success = True
    seqs = map(oagseq, R)
    stackArrays = tuple(map(partialnpchar, seqs))
    seqArray = npvstack(stackArrays)

    # print(repr(seqArray))
    quals = nparray([[chr2ph[i] for i in qual] for
                     qual in map(oagqual, R)])
    # Qualities of 2 are placeholders and mean nothing in Illumina sequencing.
    # Let's turn them into what they should be: nothing.
    # quals[quals < 3] = 0
    # --- Actually ---,  it seems that the q scores of 2 are higher quality
    # than Illumina expects them to be, so let's keep them. They should also
    # ideally be recalibrated.
    qualA = ccopy(quals)
    qualC = ccopy(quals)
    qualG = ccopy(quals)
    qualT = ccopy(quals)
    qualA[seqArray != "A"] = 0
    qualAFlat = nsum(qualA, 0, dtype=np.int64)
    qualC[seqArray != "C"] = 0
    qualCFlat = nsum(qualC, 0, dtype=np.int64)
    qualG[seqArray != "G"] = 0
    qualGFlat = nsum(qualG, 0, dtype=np.int64)
    qualT[seqArray != "T"] = 0
    qualTFlat = nsum(qualT, 0, dtype=np.int64)
    qualAllSum = npvstack(
        [qualAFlat, qualCFlat, qualGFlat, qualTFlat])
    newSeq = nparray([letterNumDict[i] for i in npargmax(qualAllSum, 0)])
    MaxPhredSum = npamax(qualAllSum, 0)  # Avoid calculating twice.
    FA = nparray([sum([seq[i] == newSeq[i] for
                       seq in seqs]) for
                  i in xrange(lenSeq)], dtype=np.int64)
    # Sums the quality score for all bases, then scales it by the number of
    # agreed bases. There could be more informative ways to do so, but
    # this is primarily a placeholder.
    phredQuals = ndiv(nmul(FA, nsum(nparray([[chr2ph[i] for i in qual] for
                                             qual in map(oagqual, R)]),
                                    0)), lenR, dtype=np.int64)
    ND = lenR * lenSeq - nsum(FA)
    # newSeq[phredQuals == 0] = "N"
    phredQuals[phredQuals < 0] = 0
    try:
        PVString = "|PV=%s" % ",".join([int2Str[i] for i in phredQuals])
    except KeyError:
        PVString = "|PV=%s" % ",".join(phredQuals.astype(str))
    try:
        phredQualsStr = "".join([ph2chrDict[i] for i in phredQuals])
    except KeyError:
        phredQualsStr = "".join(map(ph2chr, phredQuals))
    try:
        FAString = "|FA=%s" % ",".join([int2Str[i] for i in FA])
    except KeyError:
        FAString = ",".join(FA.astype(str))
    TagString = "|FM=%s|ND=%s%s%s" % (lenR, ND, FAString,
                                      PVString)
    consolidatedFqStr = "@%s %s%s\n%s\n+\n%s\n" % (name, R[0].comment,
                                                   TagString,
                                                   newSeq.tostring(),
                                                   phredQualsStr)
    if(not Success):
        return consolidatedFqStr.replace("Pass", "Fail")
    return consolidatedFqStr


@cython.returns(cystr)
def CutadaptPaired(cystr fq1, cystr fq2,
                   p3Seq="default", p5Seq="default",
                   int overlapLen=6, cython.bint makeCall=True):
    """
    Returns a string which can be called for running cutadapt v.1.7.1
    for paired-end reads in a single call.
    """
    cdef cystr outfq1
    cdef cystr outfq2
    cdef cystr commandStr
    outfq1 = ".".join(fq1.split('.')[0:-1] + ["cutadapt", "fastq"])
    outfq2 = ".".join(fq2.split('.')[0:-1] + ["cutadapt", "fastq"])
    if(p3Seq == "default"):
        raise Tim("3-prime primer sequence required for cutadapt!")
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


@cython.locals(overlapLen=int)
@cython.returns(cystr)
def CutadaptString(fq, p3Seq="default", p5Seq="default", overlapLen=6):
    """
    Returns a string which can be called for running cutadapt v.1.7.1.
    """
    outfq = ".".join(fq.split('.')[0:-1] + ["cutadapt", "fastq"])
    if(p3Seq == "default"):
        raise Tim("3-prime primer sequence required for cutadapt!")
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


@cython.locals(overlapLen=int)
def CutadaptSingle(fq, p3Seq="default", p5Seq="default", overlapLen=6):
    """
    Calls cutadapt to remove adapter sequence at either end of the reads.
    Written for v1.7.1 and single-end calls.
    """
    commandStr, outfq = CutadaptString(fq, p3Seq=p3Seq, p5Seq=p5Seq,
                                       overlapLen=overlapLen)
    subprocess.check_call(commandStr)
    return outfq


@cython.locals(overlapLap=int, numChecks=int)
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
            raise subprocess.CalledProcessError(
                fq2Popen.returncode, fq2Str, "Cutadapt failed for read 2!")


@cython.locals(useGzip=cython.bint, hpLimit=int)
def FastqPairedShading(fq1, fq2, indexFq="default",
                       useGzip=False, SetSize=10,
                       int head=2):
    """
    TODO: Unit test for this function.
    Tags fastqs with barcodes from an index fastq.
    """
    #  C declarations
    cdef pFastqProxy_t read1
    cdef pFastqProxy_t read2
    cdef pysam.cfaidx.FastqProxy indexRead
    cdef cystr outfq1
    cdef cystr outfq2
    pl("Now beginning fastq marking: Pass/Fail and Barcode")
    if(indexFq == "default"):
        raise ValueError("For an i5/i7 index ")
    outfq1 = TrimExt(fq1).replace(".fastq",
                                  "").split("/")[-1] + ".shaded.fastq"
    outfq2 = TrimExt(fq2).replace(".fastq",
                                  "").split("/")[-1] + ".shaded.fastq"
    if(useGzip):
        outfq1 += ".gz"
        outfq2 += ".gz"
    pl("Output fastqs: {}, {}.".format(outfq1, outfq2))
    inFq1 = pFastqFile(fq1)
    inFq2 = pFastqFile(fq2)
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
    inIndex = pysam.FastqFile(indexFq)
    hpLimit = len(inIndex.next().sequence) * 5 // 6
    inIndex = pysam.FastqFile(indexFq)
    ifin = inIndex.next
    outFqSet1 = []
    outFqSet2 = []
    numWritten = 0
    ofh1w = outFqHandle1.write
    ofh2w = outFqHandle2.write
    for read1 in inFq1:
        if(numWritten >= SetSize):
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
        read2 = ifn2()
        indexRead = ifin()
        tempBar = "%s%s%s" % (read1.sequence[:head], indexRead.sequence,
                              read2.sequence[:head])
        # bLen - 10 of 12 in a row, or 5/6. See Loeb, et al.
        # This is for removing low complexity reads
        # print("bLen is {}".format(bLen))
        if(BarcodePasses(tempBar, hpLimit=hpLimit)):
            tagStr = "|FP=IndexPass|BS=" + tempBar
        else:
            tagStr = "|FP=IndexFail|BS=" + tempBar
        read1.comment += tagStr
        read2.comment += tagStr
        f1.write(str(read1))
        f2.write(str(read2))
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
                       indexFq="default",
                       outfq="default",
                       cython.bint gzip=False,
                       int head=0):
    """
    Unit test done. Marks a single fastq with its index fq string.
    """
    cdef pysam.cfaidx.FastqProxy read1
    cdef pFastqProxy_t pRead1, pIndexRead
    cdef pysam.cfaidx.FastqFile inFq1, inIndex
    pl("Now beginning fastq marking: Pass/Fail and Barcode")
    if(indexFq == "default"):
        raise ValueError("For an i5/i7 index ")
    if(outfq == "default"):
        outfq = '.'.join(fq.split('.')[0:-1]) + '.shaded.fastq'
    inFq1 = pysam.FastqFile(fq)
    outFqHandle1 = open(outfq, "w")
    inIndex = pysam.FastqFile(indexFq)
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


@cython.returns(cystr)
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
        raise Tim("YOU HAVE NO CHANCE TO SURVIVE MAKE YOUR TIME")
    return tagDict


def pairedFastqConsolidate(fq1, fq2, float stringency=0.9,
                           int SetSize=100, bint parallel=True):
    """
    TODO: Unit test for this function.
    Also, it would be nice to do a groupby() that separates read 1 and read
    2 records so that it's more pythonic, but that's a hassle.
    """
    cdef cystr outFq1, outFq2
    cdef cython.int checks
    pl("Now running pairedFastqConsolidate on {} and {}.".format(fq1, fq2))
    pl("(What that really means is that I'm running "
       "singleFastqConsolidate twice")
    if(not parallel):
        outFq1 = singleFastqConsolidate(fq1, stringency=stringency,
                                        SetSize=SetSize)
        outFq2 = singleFastqConsolidate(fq2, stringency=stringency,
                                        SetSize=SetSize)
    else:
        # Make background process command string
        cStr = ("python -c 'from MawCluster.BCFastq import singleFastqConsoli"
                "date;singleFastqConsolidate(\"%s\", stringency=" % fq2 +
                "%s, SetSize=%s);import sys;sys.exit(0)'" % (stringency,
                                                             SetSize))
        # Submit
        PFC_Call = subprocess.Popen(cStr, stderr=None, shell=True,
                                    stdout=None, stdin=None, close_fds=True)
        # Run foreground job on other read
        outFq1 = singleFastqConsolidate(fq1, stringency=stringency,
                                        SetSize=SetSize)
        checks = 0
        # Have the other process wait until it's finished.
        while PFC_Call.poll() is None and checks < 3600:
            time.sleep(1)
            checks += 1
        if(PFC_Call.returncode == 0):
            return outFq1, TrimExt(fq2) + ".cons.fastq"
        elif(PFC_Call.returncode is not None):
            raise FunctionCallException(
                cStr, ("Background singleFastqConsolidate "
                       "returned non-zero exit status."),
                shell=True)
        else:
            raise FunctionCallException(
                cStr, ("Background singleFastqConsolidate took more than an"
                       "hour longer than the other read fastq. Giving up!"),
                shell=True)
    print("Consolidation a success for both!")
    return outFq1, outFq2


def singleFastqConsolidate(cystr fq, float stringency=0.9,
                           int SetSize=100,
                           cython.bint onlyNumpy=True,
                           cython.bint skipFails=False,
                           object fn=cFRF_helper,
                           object getBS=getBS,
                           object groupby=groupby):
    cdef cystr outFq, bc4fq, ffq
    cdef pFastqFile_t inFq
    cdef list StringList
    cdef int numProc
    outFq = TrimExt(fq) + ".cons.fastq"
    pl("Now running singleFastqConsolidate on {}.".format(fq))
    inFq = pFastqFile(fq)
    outputHandle = open(outFq, 'w')
    StringList = []
    workingBarcode = ""
    numProc = 0
    ohw = outputHandle.write
    sla = StringList.append
    for bc4fq, fqRecGen in groupby(inFq, key=getBS):
        ffq = fn(list(fqRecGen), bc4fq)
        sla(ffq)
        numProc += 1
        if (numProc == SetSize):
            ohw("".join(StringList))
            StringList = []
            sla = StringList.append
            numProc = 0
            continue
    ohw("".join(StringList))
    outputHandle.flush()
    inFq.close()
    outputHandle.close()
    print("Consolidation a success for inFq: %s!" % fq)
    return outFq


def TrimHomingSingle(
        fq,
        cystr homing=None,
        cystr trimfq=None,
        int bcLen=12,
        cystr trim_err=None,
        int start_trim=1):
    """
    TODO: Unit test for this function.
    """
    cdef pysam.cfaidx.FastqProxy read
    cdef int HomingLen, TotalTrim
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
        tw(SliceFastqProxy(read, firstBase=TotalTrim,
                              addString="|BS=" + read.sequence[0:bcLen]))
    trimHandle.close()
    errHandle.close()
    return trimfq


def TrimHomingPaired(inFq1, inFq2, int bcLen=12,
                     cystr homing=None, cystr trimfq1=None,
                     cystr trimfq2=None, int start_trim=1):
    """
    TODO: Unit test for this function.
    """
    cdef pysam.cfaidx.FastqProxy read1, read2
    cdef int HomingLen, TotalTrim
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
        tw1(SliceFastqProxy(read1, firstBase=TotalTrim,
                               addString="|BS=%s" % barcode))
        tw2(SliceFastqProxy(read2, firstBase=TotalTrim,
                               addString="|BS=%s" % barcode))
    trimHandle1.close()
    trimHandle2.close()
    errHandle.close()
    return trimfq1, trimfq2


@cython.locals(asDict=cython.bint)
@cython.returns(cystr)
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


@cython.returns(tuple)
def BarcodeRescueDicts(cystr indexFqPath, int minFam=10, int n=1,
                       cystr tmpFile=None):
    """Returns two dictionaries.
    1. rescueHistDict maps random barcodes to the central barcode they
    should have been.
    2. TrueFamDict's keys are all "True family" barcodes with the value set
    to None. It's just O(1) to check for a hashmap key's membership vs. O(n)
    for checking a list's membership.
    """
    cdef list histList
    cdef dict histDict, rescueHistDict, TrueFamDict
    cdef cystr h, y, x
    cdef int y1
    cStr = ("zcat %s | paste - - - - | cut -f2 | sort | " % indexFqPath +
            "uniq -c | awk 'BEGIN {{OFS=\"\t\"}};{{print $1, $2}}'")
    pl("Calling cStr: %s" % cStr)
    if tmpFile is None:
        histList = [tuple(tmpStr.split("\t")) for tmpStr in
                    check_output(cStr,
                                 shell=True).split("\n") if tmpStr != ""]
    else:
        check_call(cStr + " > %s" % tmpFile, shell=True)
        histList = [tuple(tmpStr.split("\t")) for tmpStr in
                    open(tmpFile, "r").read().split("\n") if tmpStr != ""]
    histDict = {y: int(x) for x, y in histList}
    TrueFamDict = {x: None for x, y1 in histDict.iteritems() if y1 >= minFam}
    if(len(TrueFamDict) == 0):
        rescueHistDict = {}
    else:
        rescueHistDict = {h: y for y in TrueFamDict.iterkeys() for
                          h in hamming_cousins(y, n=n)}
    if(tmpFile is not None):
        check_call(["rm", tmpFile])
    return rescueHistDict, TrueFamDict


cdef bint BarcodePasses(cystr barcode, int hpLimit=-1, bint useRe=True):
    if(hpLimit < 0):
        raise Tim("Barcode length must be set to test if it passes!")
    if(useRe):
        b = re.compile("(ACGT){%s}" % hpLimit)
        if b.match(barcode) is not None:
            return False
        if("N" in barcode):
            return False
    else:
        if("A" * hpLimit in barcode or "C" * hpLimit in barcode or
           "G" * hpLimit in barcode or "T" * hpLimit in barcode or
           "N" in barcode):
            return False
    return True


def RescueShadingWrapper(cystr inFq1, cystr inFq2, cystr indexFq=None,
                         int minFam=10, int mm=1, int head=2):
    cdef dict rescueDict, TrueFamDict
    cdef cystr tmpFilename
    tmpFilename = str(uuid.uuid4().get_hex()[0:8]) + ".tmp"
    pl("Calling RescueShadingWrapper.")
    if(indexFq is None):
        raise Tim("Index Fq must be set to rescue barcodes.")
    pl("About to do a rescue step for barcodes for %s and %s." % (inFq1,
                                                                  inFq2))
    rescueDict, TrueFamDict = BarcodeRescueDicts(indexFq, n=mm, minFam=minFam,
                                                 tmpFile=tmpFilename)
    pl("Dictionaries filled!")
    return RescuePairedFastqShading(inFq1, inFq2, indexFq,
                                    rescueDict=rescueDict,
                                    TrueFamDict=TrueFamDict, head=head)


@cython.returns(tuple)
def RescuePairedFastqShading(cystr inFq1, cystr inFq2,
                             cystr indexFq,
                             dict rescueDict=None, dict TrueFamDict=None,
                             cystr outFq1=None, cystr outFq2=None,
                             int head=2):
    """Rescues orphans from a barcode rescue.
    Works under the assumption that the number of random nucleotides used
    as molecular barcodes is sufficiently high that any read "family" with
    size below the minFam used to create the rescueDict object can be safely
    considered to be a sequencer error.
    """
    cdef cystr tmpBS, saltedBS, tagStr
    cdef pFastqFile_t inHandle1, inHandle2, indexHandle
    cdef pFastqProxy_t rec1, rec2, index_read
    cdef int bLen, hpLimit
    if(outFq1 is None):
        outFq1 = TrimExt(inFq1) + ".rescued.shaded.fastq"
    if(outFq2 is None):
        outFq2 = TrimExt(inFq2) + ".rescued.shaded.fastq"
    if(TrueFamDict is None or rescueDict is None):
        raise Tim("TrueFamDict and rescueDict must not be None! Reprs! "
                  "TFD: %s. Rescue: %s." % (repr(TrueFamDict),
                                            repr(rescueDict)))
    print("Beginning RescuePairedFastqShading. Outfqs: %s, %s." % (outFq1,
                                                                   outFq2))
    inHandle1 = pFastqFile(inFq1)
    inHandle2 = pFastqFile(inFq2)
    indexHandle = pFastqFile(indexFq)
    outHandle1 = open(outFq1, "w")
    outHandle2 = open(outFq2, "w")
    ohw1 = outHandle1.write
    ohw2 = outHandle2.write
    ih2n = inHandle2.next
    bLen = len(indexHandle.next().sequence) + 2 * head
    hpLimit = bLen * 5 // 6  # Homopolymer limit
    indexHandle = pFastqFile(indexFq)
    ihn = indexHandle.next
    for rec1 in inHandle1:
        try:
            index_read = ihn()
            rec2 = ih2n()
        except StopIteration:
            raise Tim("Index fastq and read fastqs have different sizes. "
                      "Abort!")
        tmpBS = index_read.sequence
        try:
            TrueFamDict[tmpBS]
            saltedBS = "%s%s%s" % (rec1.sequence[:head], tmpBS,
                                   rec2.sequence[:head])
            if(BarcodePasses(saltedBS, hpLimit=hpLimit)):
                tagStr = " |FP=IndexPass|BS=%s" % saltedBS
                rec1.comment += tagStr
                rec2.comment += tagStr
            else:
                tagStr = " |FP=IndexFail|BS=%s" % saltedBS
                rec1.comment += tagStr
                rec2.comment += tagStr
            ohw1(str(rec1))
            ohw2(str(rec1))
            continue
        except KeyError:
            pass
        try:
            saltedBS = rescueDict[tmpBS]
            saltedBS = rec1.sequence[:head] + saltedBS + rec2.sequence[:head]
            if(BarcodePasses(saltedBS, hpLimit=hpLimit)):
                tagStr = "|FP=IndexPass|BS=" + saltedBS + "|OS=" + tmpBS
                rec1.comment += tagStr
                rec2.comment += tagStr
            else:
                tagStr = "|FP=IndexFail|BS=" + saltedBS + "|OS=" + tmpBS
                rec1.comment += tagStr
                rec2.comment += tagStr
            ohw1(str(rec1))
            ohw2(str(rec1))
            continue
        except KeyError:
            # This isn't in a true family. Blech!
            saltedBS = rec1.sequence[:head] + tmpBS + rec2.sequence[:head]
            tagStr = "|FP=IndexFail|BS=" + saltedBS
            rec1.comment += tagStr
            rec2.comment += tagStr
            ohw1(str(rec1))
            ohw2(str(rec1))
            continue
    outHandle1.close()
    outHandle2.close()
    return outFq1, outFq2
