# cython: boundscheck=False, c_string_type=str, c_string_encoding=ascii
# cython: cdivision=True, cdivision_warnings=True, profile=True

"""
Contains various utilities for working with barcoded fastq files.
TODO: Mark family sizes not before family consolidation but during.

c_string_type and c_string_encoding should speed up the text processing.
"""

import logging
import os
import shlex
import subprocess
import copy
import gzip
import sys
import collections
import time
import cStringIO
import operator
from operator import attrgetter as oag
from operator import add as oadd
from operator import le as ole
from operator import ge as oge
from operator import div as odiv
from operator import mul as omul
from operator import add as oadd
from subprocess import check_call

import Bio
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq as BioSeq
import cython
cimport cython
import numpy as np
cimport numpy as np
from numpy import array as nparray
from numpy import less as nless
from numpy import sum as nsum
from numpy import amax as npamax
from numpy import argmax as npargmax
from numpy import multiply as nmultiply
from numpy import subtract as nsubtract
from numpy import any as npany
from numpy import vstack as npvstack
from numpy import greater as ngreater
import pysam
cimport pysam.cfaidx
from pysam.cfaidx import FastqProxy as cFastqProxy
import numconv

ctypedef np.int64_t dtypei_t

from utilBMF.HTSUtils import printlog as pl
from utilBMF.HTSUtils import PipedShellCall
from utilBMF.HTSUtils import ph2chr
from utilBMF.HTSUtils import ph2chrDict
from utilBMF.HTSUtils import chr2ph
from utilBMF import HTSUtils
from utilBMF.ErrorHandling import ThisIsMadness
from utilBMF.HTSUtils import FacePalm


letterNumDict = {}
letterNumDict['A'] = 0
letterNumDict['C'] = 1
letterNumDict['G'] = 2
letterNumDict['T'] = 3
letterNumDict[0] = 'A'
letterNumDict[1] = 'C'
letterNumDict[2] = 'G'
letterNumDict[3] = 'T'



def chr2phFunc(x):
    return chr2ph[x]


class pFastqProxy:
    """
    Python container for pysam.cfaidx.FastqProxy with persistence.
    """
    def __init__(self, FastqProxy):
        try:
            assert isinstance(FastqProxy, cFastqProxy)
        except AssertionError:
            FacePalm("pFastqProxy requires a pysam.cfaidx.FastqProxy "
                     "object for initialization!")
        self.comment = FastqProxy.comment
        self.quality = FastqProxy.quality
        self.sequence = FastqProxy.sequence
        self.name = FastqProxy.name

    def __str__(self):
        return "\n".join(["".join(["@", self.name, " ", self.comment]),
                         self.sequence,
                         "+", self.quality, ""])


@cython.locals(checks=cython.int, highMem=cython.bint,
               parallel=cython.bint)
def BarcodeSortBoth(inFq1, inFq2, sortMem="6G", parallel=True):
    if(parallel is False):
        pl("Parallel barcode sorting is set to false. Performing serially.")
        return BarcodeSort(inFq1), BarcodeSort(inFq2)
    outFq1 = '.'.join(inFq1.split('.')[0:-1] + ["BS", "fastq"])
    outFq2 = '.'.join(inFq2.split('.')[0:-1] + ["BS", "fastq"])
    pl("Sorting {} and {} by barcode sequence.".format(inFq1, inFq2))
    highMemStr = "-S " + sortMem
    if(inFq1.endswith(".gz")):
        BSstring1 = ("zcat " + inFq1 + " | paste - - - - | sed "
                     "'s: #G~:\t#G~:g' |  awk 'BEGIN {{FS=OFS=\"\t\"}};{{print"
                     " $3,$0}}' | sort -k1,1 %s | cut -f2- | " % highMemStr +
                     "sed 's:\t#G~: #G~:g' | tr '\t' '\n' > " + outFq1)
    else:
        BSstring1 = ("cat " + inFq1 + " | paste - - - -  | sed "
                     "'s: #G~:\t#G~:g' |  awk 'BEGIN {{FS=OFS=\"\t\"}};{{print"
                     " $3,$0}}' | sort -k1,1 %s | cut -f2- | " % highMemStr +
                     "sed 's:\t#G~: #G~:g' | tr '\t' '\n' > " + outFq1)
    if(inFq2.endswith(".gz")):
        BSstring2 = ("zcat " + inFq2 + " | paste - - - - | sed "
                     "'s: #G~:\t#G~:g' |  awk 'BEGIN {{FS=OFS=\"\t\"}};{{print"
                     " $3,$0}}' | sort -k1,1 %s | cut -f2- | " % highMemStr +
                     "sed 's:\t#G~: #G~:g' | tr '\t' '\n' > " + outFq2)
    else:
        BSstring2 = ("cat " + inFq2 + " | paste - - - - | sed "
                     "'s: #G~:\t#G~:g' |  awk 'BEGIN {{FS=OFS=\"\t\"}};{{print"
                     " $3,$0}}' | sort -k1,1 %s | cut -f2- | " % highMemStr +
                     "sed 's:\t#G~: #G~:g' | tr '\t' '\n' > " + outFq2)
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
                checks = oadd(checks, 1)
        raise subprocess.CalledProcessError("Barcode first sort didn't work, "
                                            "it seems - took more than 5 minu"
                                            "tes longer than the second barco"
                                            "de sort.")


@cython.locals(highMem=cython.bint)
def BarcodeSort(inFastq, outFastq="default", highMem=True):
    pl("Sorting {} by barcode sequence.".format(inFastq))
    highMemStr = ""
    if(highMem):
        highMemStr = " -S 6G "
    if(outFastq == "default"):
        outFastq = '.'.join(inFastq.split('.')[0:-1] + ["BS", "fastq"])
    if(inFastq.endswith(".gz")):
        BSstring = ("zcat " + inFastq + " | paste - - - - | sed "
                    "'s: #G~:\t#G~:g' |  awk 'BEGIN {{FS=OFS=\"\t\"}};{{print"
                    " $3,$0}}' | sort -k1,1 %s | cut -f2- | " % highMemStr +
                    "sed 's:\t#G~: #G~:g' | tr '\t' '\n' > " + outFastq)
    else:
        BSstring = ("cat " + inFastq + " | paste - - - - | sed "
                    "'s: #G~:\t#G~:g' |  awk 'BEGIN {{FS=OFS=\"\t\"}};{{print"
                    " $3,$0}}' | sort -k1,1 %s | cut -f2- | " % highMemStr +
                    "sed 's:\t#G~: #G~:g' | tr '\t' '\n' > " + outFastq)
    PipedShellCall(BSstring)
    pl("Barcode Sort shell call: {}".format(BSstring))
    # pl("Command: {}".format(BSstring.replace(
    #    "\t", "\\t").replace("\n", "\\n")))
    return outFastq


@cython.locals(stringency=cython.float, hybrid=cython.bint,
               famLimit=cython.int, keepFails=cython.bint,
               Success=cython.bint, PASS=cython.bint, frac=cython.float,
               lenR=cython.int)
def compareFastqRecords(R, stringency=0.9, famLimit=200,
                        keepFails=True):
    """
    Compares the fastq records to create a consensus sequence (if it
    passes a filter)
    """
    cdef np.ndarray[dtypei_t, ndim = 1] numAgreed
    cdef np.ndarray[dtypei_t, ndim = 1] probs
    lenR = len(R)
    try:
        famLimit = int(famLimit)
    except ValueError:
        pl("famLimit arg must be integer. Set to default: 200.")
    if(oge(lenR, famLimit)):
        logging.debug(
            "Read family - {} with {} members was capped at {}. ".format(
                R[0], len(R), famLimit))
        R = R[:famLimit]
    seqs = [str(r.seq) for r in R]
    maxScore = 0
    for seq in seqs:
        # print("Seq: {}".format(str(seq)))
        numEq = sum([seq == seqItem for seqItem in seqs])
        if(numEq > maxScore):
            maxScore = numEq
            finalSeq = seq
    frac = numEq * 1.0 / lenR
    if(ole(frac, 0.5)):
        R[0].description = R[0].description.replace("Pass", "Fail")
    else:
        return compareFastqRecordsInexactNumpy(R)
    probs = nmultiply(lenR, R[0].letter_annotations['phred_quality'],
                        dtype=np.int64)
    numAgreed = nparray([sum([seq[i] == finalSeq[i] for seq in seqs])
                          for i in range(len(finalSeq))], dtype=np.int64)
    TagString = "".join([" #G~FA=", ",".join(numAgreed.astype(str)),
                         " #G~FM=", str(lenR),
                         " #G~PV=", ",".join(probs.astype(str))])
    probs[probs > 93] = 93
    QualString = "".join([ph2chrDict[i] for i in probs])
    consFqString = "\n".join(
        ["".join(["@", R[0].description, TagString]),
         finalSeq, "+", QualString])
    return consFqString


@cython.locals(Success=cython.bint)
def compareFastqRecordsInexactNumpy(R):
    """
    Calculates the most likely nucleotide
    at each position and returns the joined record.
    """

    seqs = nparray([str(record.seq) for record in R], dtype=np.str_)
    stackArrays = tuple([np.char.array(s, itemsize=1) for s in seqs])
    seqArray = npvstack(stackArrays)
    cdef np.ndarray[dtypei_t, ndim = 2] quals
    cdef np.ndarray[dtypei_t, ndim = 2] qualA
    cdef np.ndarray[dtypei_t, ndim = 2] qualC
    cdef np.ndarray[dtypei_t, ndim = 2] qualG
    cdef np.ndarray[dtypei_t, ndim = 2] qualT
    cdef np.ndarray[dtypei_t, ndim = 1] qualAFlat
    cdef np.ndarray[dtypei_t, ndim = 1] qualCFlat
    cdef np.ndarray[dtypei_t, ndim = 1] qualGFlat
    cdef np.ndarray[dtypei_t, ndim = 1] qualTFlat
    cdef np.ndarray[dtypei_t, ndim = 2] qualAllSum
    cdef np.ndarray[dtypei_t, ndim = 1] MaxPhredSum
    cdef np.ndarray[dtypei_t, ndim = 1] phredQuals
    cdef np.ndarray[dtypei_t, ndim = 1] numFamAgreed
    # print(repr(seqArray))
    quals = nparray([
        record.letter_annotations['phred_quality']
        for record in R], dtype=np.int64)
    qualA = copy.copy(quals)
    qualC = copy.copy(quals)
    qualG = copy.copy(quals)
    qualT = copy.copy(quals)
    qualA[seqArray != "A"] = 0
    qualAFlat = nsum(qualA, 0)
    qualC[seqArray != "C"] = 0
    qualCFlat = nsum(qualC, 0)
    qualG[seqArray != "G"] = 0
    qualGFlat = nsum(qualG, 0)
    qualT[seqArray != "T"] = 0
    qualTFlat = nsum(qualT, 0)
    qualAllSum = npvstack(
        [qualAFlat, qualCFlat, qualGFlat, qualTFlat])
    finalSeq = "".join([letterNumDict[i] for i in npargmax(qualAllSum, 0)])
    MaxPhredSum = npamax(qualAllSum)  # Avoid calculating twice.
    phredQuals = nsubtract(
        nmultiply(2, MaxPhredSum),
        nsum(qualAllSum, 0))
    phredQuals[phredQuals == 0] = 93
    phredQuals[phredQuals < 0] = 0
    numFamAgreed = nparray([sum([seq[i] == finalSeq[i] for seq in seqs]) for
                             i in range(len(seqs[0]))], dtype=np.int64)
    TagString = "".join([" #G~FA=",
                         ",".join(numFamAgreed.astype(str)),
                         " #G~FM=", str(len(R)),
                         " #G~PV=", ",".join(phredQuals.astype(str))])
    phredQuals[phredQuals > 93] = 93
    QualString = "".join([ph2chrDict[i] for i in phredQuals])
    consFqString = "\n".join(
        ["".join(["@", R[0].description, TagString]), finalSeq,
         "+", QualString])
    return consFqString


@cython.locals(stringency=cython.float, hybrid=cython.bint,
               famLimit=cython.int, keepFails=cython.bint,
               Success=cython.bint, PASS=cython.bint, frac=cython.float,
               compressB64=cython.bint, lenR=cython.int,
               numEq=cython.int, maxScore=cython.int, ND=cython.int)
def compareFqRecsFqPrx(R, stringency=0.9, hybrid=False,
                       famLimit=1000, keepFails=True,
                       makeFA=True, makePV=True):
    """
    Compares the fastq records to create a consensus sequence (if it
    passes a filter)
    """
    cdef np.ndarray[dtypei_t, ndim = 1] phredQuals
    cdef np.ndarray[dtypei_t, ndim = 1] FA
    compress85 = True
    lenR = len(R)
    lenRStr = str(lenR)
    if(lenR > famLimit):
        logging.debug(
            "Read family - {} with {} members was capped at {}. ".format(
                R[0], lenR, famLimit))
        R = R[:famLimit]
    seqs = map(oag("sequence"), R)
    maxScore = 0
    Success = False
    numEq = 0
    for seq in seqs:
        numEq = sum(seq == seqItem for seqItem in seqs)
        if(oge(numEq, maxScore)):
            maxScore = numEq
            finalSeq = seq
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
        return compareFqRecsFast(R, makePV=makePV, makeFA=makeFA)
    FA = nparray([sum([seq[i] == finalSeq[i] for seq in seqs]) for i
                   in range(len(finalSeq))], dtype=np.int64)
    phredQuals = nparray([chr2ph[i] for i in list(R[0].quality)],
                          dtype=np.int64)
    phredQuals[phredQuals < 3] = 0
    phredQuals = nmultiply(lenR, phredQuals, dtype=np.int64)
    if(npany(ngreater(phredQuals, 93))):
        QualString = "".join(map(ph2chr, phredQuals))
        PVString = oadd(" #G~PV=",
                                ",".join(phredQuals.astype(str).tolist()))
    else:
        QualString = "".join([ph2chrDict[i] for i in phredQuals])
        PVString = oadd(" #G~PV=",
                                ",".join(phredQuals.astype(str).tolist()))
    TagString = "".join([" #G~FM=", lenRStr, " #G~FA=",
                         ",".join(nparray(FA).astype(str)),
                         " #G~ND=", str(nsubtract(lenR * len(seqs[0]),
                                                    nsum(FA))),
                         PVString])
    try:
        consFqString = "\n".join(
            ["".join(["@", R[0].name, " ", R[0].comment, TagString]),
             finalSeq,
             "+",
             QualString])
    except TypeError:
        print("TagString: {}".format(TagString))
        print("finalSeq: {}".format(finalSeq))
        print("QualString: {}".format(QualString))
        print("Name {}".format(R[0].name))
        print("Comment {}".format(R[0].comment))
        print("".join(["@", R[0].name, " ", R[0].comment, TagString]))
        raise ThisIsMadness("I can't figure out what's going on.")
    if(Success is False):
        return consFqString.replace("Pass", "Fail")
    return consFqString


@cython.locals(Success=cython.bint, ND=cython.int, lenR=cython.int)
def compareFqRecsFast(R, makePV=True, makeFA=True):
    """
    Calculates the most likely nucleotide
    at each position and returns the joined record string.
    """
    lenR = len(R)
    Success = True
    seqs = nparray([record.sequence for record in R], dtype=np.str_)
    stackArrays = tuple([np.char.array(s, itemsize=1) for s in seqs])
    seqArray = npvstack(stackArrays)

    cdef np.ndarray[dtypei_t, ndim = 2] quals
    cdef np.ndarray[dtypei_t, ndim = 2] qualA
    cdef np.ndarray[dtypei_t, ndim = 2] qualC
    cdef np.ndarray[dtypei_t, ndim = 2] qualG
    cdef np.ndarray[dtypei_t, ndim = 2] qualT
    cdef np.ndarray[dtypei_t, ndim = 1] qualAFlat
    cdef np.ndarray[dtypei_t, ndim = 1] qualCFlat
    cdef np.ndarray[dtypei_t, ndim = 1] qualGFlat
    cdef np.ndarray[dtypei_t, ndim = 1] qualTFlat
    cdef np.ndarray[dtypei_t, ndim = 2] qualAllSum
    cdef np.ndarray[dtypei_t, ndim = 1] MaxPhredSum
    cdef np.ndarray[dtypei_t, ndim = 1] phredQuals
    cdef np.ndarray[dtypei_t, ndim = 1] FA

    # print(repr(seqArray))
    quals = nparray(
        [map(chr2phFunc, list(record.quality)) for record in R],
        dtype=np.int64)
    # Qualities of 2 are placeholders and mean nothing in Illumina sequencing.
    # Let's turn them into what they should be: nothing.
    quals[quals < 3] = 1
    qualA = copy.copy(quals)
    qualC = copy.copy(quals)
    qualG = copy.copy(quals)
    qualT = copy.copy(quals)
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
    newSeq = "".join([letterNumDict[i] for i in npargmax(qualAllSum, 0)])
    MaxPhredSum = npamax(qualAllSum, 0)  # Avoid calculating twice.
    phredQuals = nsubtract(nmultiply(2, MaxPhredSum, dtype=np.int64),
                             nsum(qualAllSum, 0, dtype=np.int64),
                             dtype=np.int64)
    FA = nparray([sum([seq[i] == newSeq[i] for seq in seqs]) for i
                   in range(len(newSeq))], dtype=np.int64)
    ND = nsubtract(lenR * len(seqs[0]), nsum(FA))
    if(npany(nless(phredQuals, 0))):
        pl("repr of phredQuals %s" % repr(phredQuals), level=logging.DEBUG)
        phredQuals = abs(phredQuals)
    if(npany(ngreater(phredQuals, 93))):
        PVString = oadd(" #G~PV=",
                                ",".join(phredQuals.astype(str)))
        phredQuals[phredQuals > 93] = 93
        phredQualsStr = "".join([ph2chrDict[i] for i in phredQuals])
    else:
        phredQualsStr = "".join([ph2chrDict[i] for i in phredQuals])
        PVString = oadd(" #G~PV=",
                                ",".join(phredQuals.astype(str)))
    TagString = "".join([" #G~FM=", str(lenR), " #G~ND=",
                         str(ND), " #G~FA=", FA.astype(str), PVString])
    consolidatedFqStr = "\n".join([
        "".join(["@", R[0].name, " ", R[0].comment, TagString]),
        newSeq,
        "+",
        phredQualsStr])
    if(Success is False):
        return consolidatedFqStr.replace("Pass", "Fail")
    return consolidatedFqStr


@cython.locals(makeCall=cython.bint)
def CutadaptPaired(fq1, fq2, p3Seq="default", p5Seq="default",
                         overlapLen=6, makeCall=True):
    """
    Returns a string which can be called for running cutadapt v.1.7.1
    for paired-end reads in a single call.
    """
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


@cython.locals(overlapLen=cython.int)
def CallCutadapt(fq, p3Seq="default", p5Seq="default", overlapLen=6):
    """
    Calls cutadapt to remove adapter sequence at either end of the reads.
    Written for v1.7.1 and single-end calls.
    """
    commandStr, outfq = CutadaptString(fq, p3Seq=p3Seq, p5Seq=p5Seq,
                                       overlapLen=overlapLen)
    subprocess.check_call(commandStr)
    return outfq


@cython.locals(overlapLap=cython.int, numChecks=cython.int)
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


@cython.locals(useGzip=cython.bint, bLen=cython.int)
def FastqPairedShading(fq1, fq2, indexfq="default",
                       useGzip=False, readPairsPerWrite=10):
    """
    Tags fastqs with barcodes from an index fastq.
    """
    #  C declarations
    cdef pysam.cfaidx.FastqProxy read1
    cdef pysam.cfaidx.FastqProxy read2
    cdef pysam.cfaidx.FastqProxy indexRead
    pl("Now beginning fastq marking: Pass/Fail and Barcode")
    if(indexfq == "default"):
        raise ValueError("For an i5/i7 index ")
    outfq1 = ('.'.join(
        [i for i in fq1.split('.')[0:-1] if i != "fastq"] +
        ['shaded', 'fastq'])).split('/')[-1]
    outfq2 = ('.'.join(
        [i for i in fq2.split('.')[0:-1] if i != "fastq"] +
        ['shaded', 'fastq'])).split('/')[-1]
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
            if(useGzip is False):
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
        bLen = len(indexRead.sequence) * 5 / 6
        # bLen - 10 of 12 in a row, or 5/6. See Loeb, et al.
        # This is for removing low complexity reads
        # print("bLen is {}".format(bLen))
        if(("N" in tempBar or "A" * bLen in tempBar
                or "C" * bLen in tempBar or "G" * bLen in tempBar
                or "T" * bLen in tempBar)):
            '''
            pl("Failing barcode for read {} is {} ".format(indexRead,
                                                           tempBar),
               level=logging.DEBUG)
            '''
            f1.write("\n".join(["".join(["@", read1.name, " ", read1.comment,
                                " #G~FP=IndexFail #G~BS=", tempBar]),
                                read1.sequence,
                                "+", read1.quality, ""]))
            f2.write("\n".join(["".join(["@", read2.name, " ", read2.comment,
                                " #G~FP=IndexFail #G~BS=", tempBar]),
                                read2.sequence,
                                "+", read2.quality, ""]))
        else:
            f1.write("\n".join(["".join(["@", read1.name, " ", read1.comment,
                                " #G~FP=IndexPass #G~BS=", tempBar]),
                                read1.sequence,
                                "+", read1.quality, ""]))
            f2.write("\n".join(["".join(["@", read2.name, " ", read2.comment,
                                " #G~FP=IndexPass #G~BS=", tempBar]),
                                read2.sequence,
                                "+", read2.quality, ""]))
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
                       gzip=False):
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
    if(gzip):
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
            tagDict[pair[0]] = pair[1].split(' ')[0]
    except IndexError:
        pl("A value is stored with the #G~ tag which doesn't contain an =.")
        pl("tagSetEntries: {}".format(tagSetEntries))
        raise IndexError("Check that fastq description meets specifications.")
    # pl("Repr of tagDict is {}".format(tagDict))
    except TypeError:
        pl("tagSetEntries: {}".format(tagSetEntries))
        raise ThisIsMadness("YOU HAVE NO CHANCE TO SURVIVE MAKE YOUR TIME")
    return tagDict


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


@cython.locals(stringency=cython.float, readPairsPerWrite=cython.int,
               UsecProfile=cython.bint, onlyNumpy=cython.bint,
               numProc=cython.int, skipSingles=cython.bint,
               skipFails=cython.int)
def pairedFastqConsolidate(fq1, fq2, stringency=0.9,
                           readPairsPerWrite=100,
                           UsecProfile=False, onlyNumpy=True,
                           skipSingles=False, skipFails=False):
    if(UsecProfile):
        import cProfile
        import pstats
        pr = cProfile.Profile()
        pr.enable()
    outFqPair1 = '.'.join(fq1.split('.')[0:-1] + ["cons", "fastq"])
    outFqPair2 = '.'.join(fq2.split('.')[0:-1] + ["cons", "fastq"])
    pl("Now running pairedFastqConsolidate on {} and {}.".format(fq1, fq2))
    pl("Command required to duplicate this action:"
       " pairedFastqConsolidate('{}', '{}', ".format(fq1, fq2) +
       "stringency={}, readPairsPerWrite={})".format(stringency,
                                                     readPairsPerWrite))
    inFq1 = SeqIO.parse(fq1, "fastq")
    inFq2 = SeqIO.parse(fq2, "fastq")
    outputHandle1 = open(outFqPair1, 'w')
    outputHandle2 = open(outFqPair2, 'w')
    # cString1 = cStringIO.StringIO()
    # cString2 = cStringIO.StringIO()
    StringList1 = []
    StringList2 = []
    workingBarcode = ""
    previousBarcode = ""
    workingSet1 = []
    workingSet2 = []
    numProc = 0
    ws1a = workingSet1.append
    ws2a = workingSet2.append
    oh1w = outputHandle1.write
    oh2w = outputHandle2.write
    sl1a = StringList1.append
    sl2a = StringList2.append
    while True:
        if(numProc % readPairsPerWrite == 0):
            # outputHandle1.write(cString1.getvalue())
            # outputHandle2.write(cString2.getvalue())
            oh1w("".join(StringList1))
            oh2w("".join(StringList2))
            StringList1 = []
            StringList2 = []
            sl1a = StringList1.append
            sl2a = StringList2.append
            # cString1 = cStringIO.StringIO()
            # cString2 = cStringIO.StringIO()
        try:
            fqRec = next(inFq1)
        except StopIteration:
            break
        bc4fq1 = GetDescTagValue(fqRec.description, "BS")
        fqRec2 = next(inFq2)
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
            ws1a(fqRec)
            ws2a(fqRec2)
            continue
        elif(workingBarcode != bc4fq1):
            if(skipSingles and len(workingSet1) == 1):
                workingBarcode = ""
                workingSet1 = []
                workingSet2 = []
                ws1a = workingSet1.append
                ws2a = workingSet2.append
                continue
            # cString1.write(compareFqRecsFqPrx(workingSet1) + "\n")
            # cString2.write(compareFqRecsFqPrx(workingSet2) + "\n")
            # String1 += compareFqRecsFqPrx(workingSet1) + "\n"
            # String2 += compareFqRecsFqPrx(workingSet2) + "\n"
            tStr1 = oadd(compareFastqRecords(workingSet1), "\n")
            tStr2 = oadd(compareFastqRecords(workingSet2), "\n")
            if(skipFails and ("Fail" in tStr1 or "Fail" in tStr2)):
                continue
            sl1a(tStr1)
            sl2a(tStr2)
            workingSet1 = [fqRec]
            workingSet2 = [fqRec2]
            workingBarcode = bc4fq1
            numProc = oadd(numProc, 1)
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
    oh1w("".join(StringList1))
    oh2w("".join(StringList2))
    outputHandle1.flush()
    outputHandle2.flush()
    inFq1.close()
    inFq2.close()
    # cString1.close()
    # cString2.close()
    outputHandle1.close()
    outputHandle2.close()
    return outFqPair1, outFqPair2


@cython.locals(stringency=cython.float, readPairsPerWrite=cython.int,
               UsecProfile=cython.bint, onlyNumpy=cython.bint)
def pairedFastqConsolidateFaster(fq1, fq2, stringency=0.9,
                                 readPairsPerWrite=100, UsecProfile=False,
                                 onlyNumpy=False, skipSingles=False,
                                 skipFails=False):
    if(UsecProfile):
        import cProfile
        import pstats
        pr = cProfile.Profile()
        pr.enable()
    outFqPair1 = '.'.join(fq1.split('.')[0:-1] + ["cons", "fastq"])
    outFqPair2 = '.'.join(fq2.split('.')[0:-1] + ["cons", "fastq"])
    pl("Now running pairedFastqConsolidateFaster on {} and {}.".format(fq1,
                                                                       fq2))
    pl("Command required to duplicate this action:"
       " pairedFastqConsolidateFaster('{}', '{}', ".format(fq1, fq2) +
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
            tStr1 = oadd(compareFqRecsFqPrx(workingSet1), "\n")
            tStr2 = oadd(compareFqRecsFqPrx(workingSet2), "\n")
            if(skipFails and ("Fail" in tStr1 or "Fail" in tStr2)):
                continue
            sl1a(tStr1)
            sl2a(tStr2)
            workingSet1 = [fqRec]
            workingSet2 = [fqRec2]
            workingBarcode = bc4fq1
            numProc = oadd(numProc, 1)
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


def renameReads(fq1, fq2, outfq1="default", outfq2="default"):
    """
    Requires barcode-sorted, consolidated fastq files,
    filtered to only the shared barcode families
    """
    if(outfq1 == "default"):
        outfq1 = fq1.split('.')[0] + '.cons.R1.fastq'
    if(outfq2 == "default"):
        outfq2 = fq2.split('.')[0] + '.cons.R2.fastq'
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
        if(("N" in bc4fq or "A" * bLen in bc4fq
                or "C" * bLen in bc4fq or "G" * bLen in bc4fq or "T" * bLen)):
            continue
        if(ole(int(GetDescTagValue(fqRec.description, "FM")), 2)):
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
            BioSeq(str(rec.seq)[0:bar_len], "fastq"),
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
        post_tag = SeqRecord(BioSeq(str(rec.seq)[TotalTrim:], "fastq"),
                             id=rec.id,
                             description=descString)
        post_tag.letter_annotations['phred_quality'] = rec.letter_annotations[
            'phred_quality'][TotalTrim:]
        SeqIO.write(post_tag, trimOpen, "fastq")
    tagsOpen.close()
    trimOpen.close()
    return(tags_file, trimfq)


def compareConsSpeed(fq1, fq2, stringency=0.666, readPairsPerWrite=100):
    """
    Made this quick script to compare speed for consolidation.
    """
    import subprocess
    import time
    from time import gmtime, strftime
    print(strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime()))
    b = subprocess.Popen("cp {} CyCons.R1.fastq".format(fq1), shell=True)
    subprocess.check_call("cp {} CyCons.R2.fastq".format(fq2), shell=True)
    while True:
        if(b.poll() != 0):
            time.sleep(1)
            continue
        else:
            break
    out1, out2 = pairedFastqConsolidate(
        "CyCons.R1.fastq", "CyCons.R2.fastq", stringency=stringency,
        readPairsPerWrite=readPairsPerWrite)
    print("Time after CyCons: {}".format(
        strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())))
    b = subprocess.Popen("cp {} PyCons.R1.fastq".format(fq1), shell=True)
    subprocess.check_call("cp {} PyCons.R2.fastq".format(fq2), shell=True)
    out1, out2 = pairedFastqConsolidate(
        "CyCons.R1.fastq", "CyCons.R2.fastq", stringency=stringency,
        readPairsPerWrite=readPairsPerWrite, onlyNumpy=True)
    print("Time after CyCons (only Numpy): {}".format(
        strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())))
    return


def CalcFamUtils(inFq):
    """
    Uses bioawk to get summary statistics on a fastq file blazingly quickly.
    """
    bioawk -c fastx '{realFams += 1; readFmSum += array[2]}};END {print "MFS="sum /len";NumRealFams="realFams";ReadFamSum="readFmSum/realFams";NumSingletons="singletons}' cfDNA06-6862_A1_ACAGTG_L001_R1_001.Tiny.shaded.BS.cons.fastq
    commandStr = ("bioawk -c fastx {{n=split($comment,array," ");k=array[4];n="
                  "split(k,array,\"=\"); len += 1; sum += array[2]; if(array[2"
                  "] == 1) {{singletons += 1}};if(array[2] >= 2) {{realFams +="
                  " 1; readFmSum += array[2]}}}};END {{print \"MFS=\"sum /len"
                  "\";NumRealFams=\"realFams\";ReadFamSum=\"readFmSum/realFams"
                  "\";NumSingletons=\"singletons}}' %s" % inFq)