# cython: c_string_type=str, c_string_encoding=ascii
# cython: cdivision=True, profile=True

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
from array import array
from string import maketrans
from numpy.core.multiarray import int_asbuffer
from cPickle import load
import ctypes

import cython
import numpy as np
import pysam
from itertools import groupby
from numpy import sum as nsum

from utilBMF.HTSUtils import (SliceFastqProxy,
                              printlog as pl,
                              pFastqProxy, TrimExt, pFastqFile, getBS,
                              hamming_cousins, GetParallelDMPPopen)
from utilBMF import HTSUtils
from utilBMF.ErrorHandling import ThisIsMadness as Tim, FunctionCallException
from utilBMF.ErrorHandling import UnsetRequiredParameter, ImproperArgumentError
try:
    from re2 import compile as regex_compile
except ImportError:
    pl("Note: re2 import failed. Fell back to re.", level=logging.DEBUG)
    from re import compile as regex_compile


ARGMAX_TRANSLATE_STRING = maketrans('\x00\x01\x02\x03', 'ACGT')
nucs = array('B', [65, 67, 71, 83])
DEF size32 = 4
DEF sizedouble = 8


def SortAndMarkFastqsCommand(Fq1, Fq2, IndexFq):
    return ("pr -mts <(cat %s | paste - - - -) <(cat %s | " % (Fq1, Fq2) +
            "paste - - - -) <(cat %s | paste - - - -) | awk" % IndexFq +
            "'BEGIN {{FS=OFS=\"\t\"}};{{BS=\"|BS=\"$10\"\"substr($2,"
            "0,2)\"\"substr($6,0,2); print $1\"\"BS, $2, $3, $4, $5\"\"BS, "
            "$6, $7, $8}}' | sort -t\"|\" -k2,2 -k1,1 | tr '\t' '\n'")


@cython.locals(checks=int,
               parallel=cython.bint, sortMem=cystr)
def BarcodeSortBoth(cystr inFq1, cystr inFq2,
                    cystr sortMem=None, cython.bint parallel=False):
    cdef cystr outFq1, outFq2, highMemStr
    if(sortMem is None):
        sortMem = "6G"
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
                cystr mem="6G", int threads=1):
    cdef cystr BSstring
    pl("Sorting {} by barcode sequence.".format(inFastq))
    BSstring = getBarcodeSortStr(inFastq, outFastq=outFastq, mem=mem,
                                 threads=threads)
    check_call(BSstring, shell=True)
    pl("Barcode Sort shell call: {}".format(BSstring))
    if(outFastq == "default"):  # Added for compatibility with getBSstr
        outFastq = '.'.join(inFastq.split('.')[0:-1] + ["BS", "fastq"])
    return outFastq


@cython.returns(cystr)
def getBarcodeSortStr(inFastq, outFastq="default", mem="",
                      int threads=1):
    if(mem != ""):
        mem = " -S " + mem
    threadStr = ""
    if(threads != 1):
        import warnings
        warnings.warn("Note: getBarcodeSortStr is being given a threads argum"
                      "ent. This is deprecated for uniformity between linux "
                      "distributions. Parameter simply being ignored.",
                      DeprecationWarning)
    if(outFastq == "default"):
        outFastq = '.'.join(inFastq.split('.')[0:-1] + ["BS", "fastq"])
    if(inFastq.endswith(".gz")):
        return ("zcat %s | paste - - - - | sort -t'|' -k3,3 -k1,1" % inFastq +
                " %s %s| tr '\t' '\n' > %s" % (mem, threadStr, outFastq))
    else:
        return ("cat %s | paste - - - - | sort -t'|' -k3,3 -k1,1" % inFastq +
                " %s %s | tr '\t' '\n' > %s" % (mem, threadStr, outFastq))


cpdef cystr QualArr2QualStr(ndarray[int32_t, ndim=1] qualArr):
    """
    cpdef wrapper for QualArr2QualStr

    Compared speed for a typed list comprehension against that for a map.
    Results:
    In [6]: %timeit QualArr2PVString(c)
    100000 loops, best of 3: 6.49 us per loop

    In [7]: %timeit QualArr2PVStringMap(c)
    100000 loops, best of 3: 8.9 us per loop

    """
    return cQualArr2QualStr(qualArr)


cpdef cystr QualArr2PVString(ndarray[int32_t, ndim=1] qualArr):
    return cQualArr2PVString(qualArr)


cpdef cystr pCompareFqRecsFast(list R, cystr name=None):
    return cCompareFqRecsFast(R, name)


cdef inline ndarray[int8_t, ndim=2] cRecListTo2DCharArray(list R):
    cdef pFastqProxy_t rec
    return np.array([cs_to_ia(rec.sequence) for rec in R],
                    dtype=np.int8)


cpdef ndarray[char, ndim=2] RecListTo2DCharArray(list R):
    return cRecListTo2DCharArray(R)


cdef class Qual2DArray:
    def __cinit__(self, size_t nRecs, size_t rLen):
        self.nRecs = nRecs
        self.rLen = rLen
        self.qualities = <int32_t *>malloc(nRecs * rLen * size32)

    def __dealloc__(self):
        free(self.qualities)


cdef class SeqQual:
    def __cinit__(self, size_t length):
        self.Seq = <int8_t *>malloc(length)
        self.Agree = <int32_t *>malloc(length * size32)
        self.Qual = <int32_t *>malloc(length * size32)
        self.length = length

    def __dealloc__(self):
        free(self.Seq)
        free(self.Agree)
        free(self.Qual)


cdef class SumArraySet:
    def __cinit__(self, size_t length):
        self.length = length
        self.counts = <int32_t *>calloc(length * 4, size32)
        self.chiSums = <double_t *>calloc(length * 4, sizedouble)
        self.argmax_arr = <int8_t *>malloc(length)

    def __dealloc__(self):
        free(self.counts)
        free(self.chiSums)
        free(self.argmax_arr)


cdef SeqQual_t cFisherFlatten(int8_t * Seqs, int32_t * Quals,
                              size_t rLen, size_t nRecs):
    cdef SeqQual_t ret
    cdef SumArraySet_t Sums
    ret = SeqQual(rLen)
    Sums = SumArraySet(rLen)
    cdef py_array Seq, Qual, Agree
    # Temporary data structures
    cdef size_t query_index = 0
    cdef size_t chisum_index = 0
    cdef size_t ndIndex = 0
    cdef size_t numbases = rLen * nRecs
    cdef size_t rLen2 = 2 * rLen
    cdef size_t rLen3 = 3 * rLen
    cdef size_t offset = 0
    cdef double_t invchi
    cdef int32_t tmpQual
    while query_index < numbases:
        offset = query_index % rLen
        if(Seqs[query_index] == 67):
            # Add in the chiSum for the observed base
            ndIndex = offset + rLen
            Sums.chiSums[ndIndex] += CHI2_FROM_PHRED(Quals[query_index])
            # Add in thi CHI2 sum contribution for the bases that
            # weren't observed.
            invchi = INV_CHI2_FROM_PHRED(Quals[query_index])
            Sums.chiSums[offset] += invchi
            Sums.chiSums[offset + rLen2] += invchi
            Sums.chiSums[offset + rLen3] += invchi
            # case "C"
        elif(Seqs[query_index] == 71):
            # case "G"
            ndIndex = offset + rLen2
            Sums.chiSums[ndIndex] += CHI2_FROM_PHRED(Quals[query_index])
            invchi = INV_CHI2_FROM_PHRED(Quals[query_index])
            Sums.chiSums[offset] += invchi
            Sums.chiSums[offset + rLen] += invchi
            Sums.chiSums[offset + rLen3] += invchi
        elif(Seqs[query_index] == 84):
            # case "T"
            ndIndex = offset + rLen3
            Sums.chiSums[ndIndex] += CHI2_FROM_PHRED(Quals[query_index])
            invchi = INV_CHI2_FROM_PHRED(Quals[query_index])
            Sums.chiSums[offset] += invchi
            Sums.chiSums[offset + rLen] += invchi
            Sums.chiSums[offset + rLen2] += invchi
        elif(Seqs[query_index] == 65):
            # case "A"
            ndIndex = offset
            Sums.chiSums[offset] += CHI2_FROM_PHRED(Quals[query_index])
            invchi = INV_CHI2_FROM_PHRED(Quals[query_index])
            Sums.chiSums[offset + rLen] += invchi
            Sums.chiSums[offset + rLen2] += invchi
            Sums.chiSums[offset + rLen3] += invchi
        else:
            # case "N"
            pass
        Sums.counts[ndIndex] += 1
        query_index += 1
    query_index = 0
    while query_index < rLen:
        # Find the most probable base
        Sums.argmax_arr[query_index] = argmax4(
            Sums.chiSums[query_index],
            Sums.chiSums[query_index + rLen],
            Sums.chiSums[query_index + rLen2],
            Sums.chiSums[query_index + rLen3])
        # Convert the argmaxes to letters
        ret.Seq[query_index] = ARGMAX_CONV(Sums.argmax_arr[query_index])
        ndIndex = query_index + Sums.argmax_arr[query_index] * rLen
        # Round
        tmpQual =  <int32_t> (- 10 * c_log10(
            igamc_pvalues(Sums.counts[ndIndex],
                          Sums.chiSums[ndIndex])) + 0.5)
        # Eliminate underflow p values
        ret.Qual[query_index] = tmpQual
        # Count agreement
        ret.Agree[query_index] = Sums.counts[ndIndex]
        query_index += 1
    return ret


cpdef SeqQual_t FisherFlatten(
        ndarray[int8_t, ndim=2, mode="c"] Seqs,
        ndarray[int32_t, ndim=2, mode="c"] Quals,
        size_t rLen, size_t nRecs):
    """
    :param Quals: [ndarray[int32_t, ndim=2]/arg] - numpy 2D-array for
    holding the qualities for a set of reads.
    :param Seqs: [ndarray[int32_t, ndim=2]/arg] - numpy 2D-array for
    holding the sequences for a set of reads.
    :param rLen: [size_t/arg] Read length
    :param nRecs: [size_t/arg] Number of records to flatten
    :return: [SeqQual_t] An array of sequence and an array of qualities
    after flattening.
    """
    return cFisherFlatten(&Seqs[0,0], &Quals[0,0], rLen, nRecs)


'''

This method works (change the second argument on amax/argmax calls
on the output), but it looks like it's slower.
I might have to do this manually.

@cython.boundscheck(False)
@cython.wraparound(False)
cdef ndarray[int32_t, ndim=2] FlattenSeqs(ndarray[char, ndim=2] seqs,
                                             ndarray[int32_t, ndim=2] quals,
                                             size_t lenR):
    cdef size_t readlen = len(seqs[0])
    cdef size_t read_index, base_index
    cdef ndarray[int32_t, ndim=2] QualSumAll
    cdef int32_t tmpInt
    QualSumAll = np.zeros([readlen, 4], dtype=np.int32)
    for read_index in range(lenR):
        for base_index in range(readlen):
            tmpInt = seqs[read_index][base_index]
            if(tmpInt == 84):
                QualSumAll[base_index][3] += quals[read_index][base_index]
            elif(tmpInt == 67):
                QualSumAll[base_index][1] += quals[read_index][base_index]
            elif(tmpInt == 71):
                QualSumAll[base_index][2] += quals[read_index][base_index]
            else:
                QualSumAll[base_index][0] += quals[read_index][base_index]
    return QualSumAll
'''


@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline int32_t RESCALE_QUALITY(
        int8_t qual, int32_t * RescalingArray) nogil:
    return RescalingArray[qual - 2]


cdef ndarray[int32_t, ndim=2, mode="c"] ParseAndRescaleQuals(
        list R, size_t rLen, size_t nRecs, py_array RescalingArray):
    cdef pFastqProxy_t rec
    cdef py_array qual_fetcher
    cdef ndarray[int32_t, ndim=2, mode="c"] ret
    cdef Qual2DArray_t Rescaler
    cdef size_t offset = 0
    cdef size_t i
    Rescaler = Qual2DArray(nRecs, rLen)
    ret = np.zeros([nRecs, rLen], dtype=np.int32)
    for rec in R:
        qual_fetcher = cs_to_ph(rec.quality)
        for i in range(rLen):
            Rescaler.qualities[i + offset] = RESCALE_QUALITY(
                qual_fetcher[i], <int32_t *> RescalingArray.data.as_ints)
        offset += rLen
    memcpy(&ret[0,0], &Rescaler.qualities, size32 * rLen * nRecs)
    return ret


cpdef cystr FastFisherFlattening(list R, cystr name=None,
                                 bint cap_quality=False,
                                 int8_t cap=38):
    return cFastFisherFlattening(R, name=name, cap_quality=cap_quality,
                                 cap=cap)


@cython.boundscheck(False)
@cython.initializedcheck(False)
@cython.wraparound(False)
cdef cystr cFastFisherFlattenRescale(list R, py_array RescalingArray,
                                     cystr name=None):
    """
    TODO: Unit test for this function.
    Calculates the most likely nucleotide
    at each position and returns the joined record string.
    After inlining:
    In [21]: %timeit pCompareFqRecsFast(fam)
   1000 loops, best of 3: 518 us per loop

    In [22]: %timeit cFRF_helper(fam)
    1000 loops, best of 3: 947 us per loop

    """
    cdef int lenR, ND, lenSeq
    cdef double_t tmpFlt
    cdef cystr PVString, TagString, newSeq, phredQualsStr, FAString
    cdef cystr consolidatedFqStr
    cdef ndarray[int32_t, ndim=2, mode="c"] quals
    cdef ndarray[int8_t, ndim=2, mode="c"] seqArray
    cdef py_array Seq, Qual, Agree
    cdef pFastqProxy_t rec
    cdef SeqQual_t ret
    if(name is None):
        name = R[0].name
    lenR = len(R)
    lenSeq = len(R[0].sequence)
    if lenR == 1:
        rec = R[0]
        TagString = ("|FM=1|ND=0|FA=" + "1" + "".join([",1"] * (lenSeq - 1)) +
                     PyArr2PVString(rec.getQualArray()))
        return "@%s %s%s\n%s\n+\n%s\n" % (name, rec.comment,
                                          TagString, rec.sequence,
                                          rec.quality)
    seqArray = cRecListTo2DCharArray(R)

    quals = ParseAndRescaleQuals(R, lenSeq, lenR, RescalingArray)
    '''
    quals = np.array([cs_to_ph(rec.quality) for
                      rec in R], dtype=np.int32)
    '''
    # TODO: Speed up this copying by switching to memcpy
    # Flatten with Fisher.
    ret = FisherFlatten(seqArray, quals, lenSeq, lenR)
    # Copy out results to python arrays
    # Seq
    Seq = array('B')
    c_array.resize(Seq, lenSeq)
    memcpy(Seq.data.as_voidptr, <int8_t*> ret.Seq, lenSeq)

    # Qual
    Qual = array('i')
    c_array.resize(Qual, lenSeq)
    memcpy(Qual.data.as_voidptr, <int32_t*> ret.Qual, lenSeq * size32)

    # Agree
    Agree = array('i')
    c_array.resize(Agree, lenSeq)
    memcpy(Agree.data.as_voidptr, <int32_t*> ret.Agree, lenSeq * size32)
    # Get seq string

    newSeq = Seq.tostring()
    # Get quality strings (both for PV tag and quality field)
    phredQualsStr = PyArr2QualStr(Qual)
    PVString = PyArr2PVString(Qual)

    # Use the number agreed
    FAString = PyArr2FAString(Agree)
    ND = lenR * lenSeq - nsum(Agree)
    TagString = "|FM=%s|ND=%s" % (lenR, ND) + PVString + FAString
    consolidatedFqStr = ("@" + name + " " + R[0].comment + TagString + "\n" +
                         newSeq + "\n+\n%s\n" % phredQualsStr)
    return consolidatedFqStr


@cython.boundscheck(False)
@cython.initializedcheck(False)
@cython.wraparound(False)
cdef cystr cFastFisherFlattening(list R,
                                 cystr name=None,
                                 bint cap_quality=False,
                                 int8_t cap=38):
    """
    TODO: Unit test for this function.
    Calculates the most likely nucleotide
    at each position and returns the joined record string.
    After inlining:
    In [21]: %timeit pCompareFqRecsFast(fam)
   1000 loops, best of 3: 518 us per loop

    In [22]: %timeit cFRF_helper(fam)
    1000 loops, best of 3: 947 us per loop

    """
    cdef int lenR, ND, lenSeq
    cdef double_t tmpFlt
    cdef cystr PVString, TagString, newSeq, phredQualsStr, FAString
    cdef cystr consolidatedFqStr
    cdef ndarray[int32_t, ndim=2, mode="c"] quals
    cdef ndarray[int8_t, ndim=2, mode="c"] seqArray
    cdef py_array Seq, Qual, Agree
    cdef pFastqProxy_t rec
    cdef SeqQual_t ret
    if(name is None):
        name = R[0].name
    lenR = len(R)
    lenSeq = len(R[0].sequence)
    if lenR == 1:
        rec = R[0]
        TagString = ("|FM=1|ND=0|FA=" + "1" + "".join([",1"] * (lenSeq - 1)) +
                     PyArr2PVString(rec.getQualArray()))
        return "@%s %s%s\n%s\n+\n%s\n" % (name, rec.comment,
                                          TagString, rec.sequence,
                                          rec.quality)
    seqArray = cRecListTo2DCharArray(R)

    '''
    quals = ParseAndRescaleQuals(R)
    '''
    quals = np.array([cs_to_ph(rec.quality) for
                      rec in R], dtype=np.int32)
    # TODO: Speed up this copying by switching to memcpy
    # Flatten with Fisher.
    ret = FisherFlatten(seqArray, quals, lenSeq, lenR)
    # Copy out results to python arrays
    # Seq
    Seq = array('B')
    c_array.resize(Seq, lenSeq)
    memcpy(Seq.data.as_voidptr, <int8_t*> ret.Seq, lenSeq)

    # Qual
    Qual = array('i')
    c_array.resize(Qual, lenSeq)
    memcpy(Qual.data.as_voidptr, <int32_t*> ret.Qual, lenSeq * size32)

    # Agree
    Agree = array('i')
    c_array.resize(Agree, lenSeq)
    memcpy(Agree.data.as_voidptr, <int32_t*> ret.Agree, lenSeq * size32)
    # Get seq string
    newSeq = Seq.tostring()
    # Get quality strings (both for PV tag and quality field)
    '''
    # For debugging.
    '''
    cdef int32_t tmp
    for tmp in Qual:
        if(tmp < 0):
            sys.stderr.write("Repr of Seq: %s\n" % Seq)
            sys.stderr.write("Str of a read: %s\n" % str(R[0]))
            sys.stderr.write("Repr of Agree: %s\n" % Agree)
            sys.stderr.write("Repr of Qual: %s\n" % Qual)
    '''
    '''
    phredQualsStr = PyArr2QualStr(Qual) if(
        not cap_quality) else CappedPyArr2QualStr(Qual, cap=cap)
    PVString = PyArr2PVString(Qual)
    # Use the number agreed
    FAString = PyArr2FAString(Agree)
    ND = lenR * lenSeq - nsum(Agree)
    TagString = "|FM=%s|ND=%s" % (lenR, ND) + PVString + FAString
    consolidatedFqStr = ("@" + name + " " + R[0].comment + TagString + "\n" +
                         newSeq + "\n+\n%s\n" % phredQualsStr)
    return consolidatedFqStr


@cython.boundscheck(False)
@cython.wraparound(False)
cdef cystr cCompareFqRecsFast(list R,
                              cystr name=None,
                              double_t minPVFrac=0.1,
                              double_t minFAFrac=0.2,
                              double_t minMaxFA=0.9):
    """
    TODO: Unit test for this function.
    Calculates the most likely nucleotide
    at each position and returns the joined record string.
    After inlining:
    In [21]: %timeit pCompareFqRecsFast(fam)
   1000 loops, best of 3: 518 us per loop

    In [22]: %timeit cFRF_helper(fam)
    1000 loops, best of 3: 947 us per loop

    """
    cdef int lenR, ND, lenSeq, tmpInt, i
    cdef double_t tmpFlt
    cdef cython.bint Success
    cdef cystr PVString, TagString, newSeq
    cdef cystr consolidatedFqStr
    # cdef char tmpChar
    cdef ndarray[int32_t, ndim=2] quals, qualA, qualC, qualG
    cdef ndarray[int32_t, ndim=2] qualT, qualAllSum
    cdef ndarray[int32_t, ndim=1] qualAFlat, qualCFlat, qualGFlat, FA
    cdef ndarray[int32_t, ndim=1] phredQuals, qualTFlat
    cdef ndarray[char, ndim=2] seqArray
    cdef py_array tmpArr
    cdef pFastqProxy_t rec
    cdef char tmpChar
    if(name is None):
        name = R[0].name
    lenR = len(R)
    lenSeq = len(R[0].sequence)
    if lenR == 1:
        phredQuals = np.array(R[0].getQualArray(), dtype=np.int32)
        TagString = ("|FM=1|ND=0|FA=" + "1" + "".join([",1"] * (lenSeq - 1)) +
                     cQualArr2PVString(phredQuals))
        return "@%s %s%s\n%s\n+\n%s\n" % (name, R[0].comment,
                                          TagString, R[0].sequence,
                                          R[0].quality)
    Success = True
    '''
    stackArrays = tuple([np.char.array(rec.sequence, itemsize=1) for rec in R])
    seqArray = np.vstack(stackArrays)
    '''
    seqArray = cRecListTo2DCharArray(R)

    quals = np.array([cs_to_ph(rec.quality) for
                      rec in R], dtype=np.int32)
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
    qualA[seqArray != 65] = 0
    qualAFlat = nsum(qualA, 0, dtype=np.int32)
    qualC[seqArray != 67] = 0
    qualCFlat = nsum(qualC, 0, dtype=np.int32)
    qualG[seqArray != 71] = 0
    qualGFlat = nsum(qualG, 0, dtype=np.int32)
    qualT[seqArray != 84] = 0
    qualTFlat = nsum(qualT, 0, dtype=np.int32)
    qualAllSum = np.vstack(
        [qualAFlat, qualCFlat, qualGFlat, qualTFlat])
    tmpArr = array('B', np.argmax(qualAllSum, 0))
    newSeq = tmpArr.tostring().translate(ARGMAX_TRANSLATE_STRING)
    phredQuals = np.amax(qualAllSum, 0)  # Avoid calculating twice.

    # Filtering
    # First, kick out bad/discordant bases
    tmpFlt = float(np.max(phredQuals))
    PVFracDivisor = minPVFrac * tmpFlt
    phredQuals[phredQuals < PVFracDivisor] = 0
    # Second, flag families which are probably not really "families"
    FA = np.array([sum([rec.sequence[i] == newSeq[i] for
                        rec in R]) for
                  i in xrange(lenSeq)], dtype=np.int32)
    if(np.min(FA) < minFAFrac * lenR or np.max(FA) < minMaxFA * lenR):
        #  If there's any base on which a family agrees less often
        #  than minFAFrac, nix the whole family.
        #  Along with that, require that at least one of the bases agree
        #  to some fraction. I've chosen 0.9 to get rid of junk families.
        phredQuals[:] = 0
    # Sums the quality score for all bases, then scales it by the number of
    # agreed bases. There could be more informative ways to do so, but
    # this is primarily a placeholder.
    ND = lenR * lenSeq - nsum(FA)
    phredQuals[phredQuals < 0] = 0
    PVString = cQualArr2PVString(phredQuals)
    phredQualsStr = cQualArr2QualStr(phredQuals)
    FAString = cQualArr2FAString(FA)
    TagString = "|FM=%s|ND=%s" % (lenR, ND) + FAString + PVString
    '''
    consolidatedFqStr = "@%s %s%s\n%s\n+\n%s\n" % (name, R[0].comment,
                                                   TagString,
                                                   newSeq,
                                                   phredQualsStr)
    In [67]: %timeit omgzwtf = "@%s %s%s\n%s\n+\n%s\n" % (name, b.comment,
                                                          TagString,
                                                          newSeq,
                                                          phredQualsStr)
    1000000 loops, best of 3: 585 ns per loop
    In [68]: %timeit omgzwtf = ("@" + name + " " + b.comment + TagString +
                                "\n" + newSeq + "\n+\n%s\n" % phredQualsStr)
    1000000 loops, best of 3: 512 ns per loop
    '''
    consolidatedFqStr = ("@" + name + " " + R[0].comment + TagString + "\n" +
                         newSeq + "\n+\n%s\n" % phredQualsStr)
    if(not Success):
        return consolidatedFqStr.replace("Pass", "Fail")
    return consolidatedFqStr


@cython.returns(cystr)
def CutadaptPaired(cystr fq1, cystr fq2,
                   p3Seq="default", p5Seq="default",
                   int overlapLen=6, cython.bint makeCall=True,
                   outfq1=None, outfq2=None):
    """
    Returns a string which can be called for running cutadapt v.1.7.1
    for paired-end reads in a single call.
    """
    cdef cystr commandStr
    if(outfq1 is None):
        outfq1 = ".".join(fq1.split('.')[0:-1] + ["cutadapt", "fastq"])
    if(outfq2 is None):
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


def DispatchParallelDMP(fq1, fq2, indexFq="default",
                        int head=-1, sortMem=None, overlapLen=None,
                        int threads=-1, int num_nucs=-1,
                        cystr p3Seq=None, cystr p5Seq=None,
                        bint cap_quality=False, int8_t cap=38):
    """
    DispatchParallelDMP prepares a PopenDispatcher instance for demultiplexing
    barcoded molecular families in parallel.
    :param fq1 [cystr/arg] - Path to read 1 fastq
    :param fq2 [cystr/arg] - Path to read 2 fastq
    :param indexFq [cystr/kwarg/"default"] - Path to index fastq
    :param head [int/kwarg/-1] - Number of bases from each of
    read 1 and read 2 with which to "salt" the barcode. Required.
    :param sortMem [cystr/kwarg/None] - sort memory argument for each
    slave bmftools dmp process. Defaults to 768M
    :param overlapLen [object/kwarg/None] - If SlaveDMPCommandString receives
    a None value for overlapLen, it defaults to 6 for the -O option for
    cutadapt.
    :param threads [int/kwarg/-1] - Required to be > 0
    :param num_nucs [int/kwarg/-1] - Required to be > 0. Number of initial
    bases to use to split reads into fastq files by barcode start.
    :param p3Seq [cystr/kwarg/None] - 3' adapter sequence
    :param p5Seq [cystr/kwarg/None] - 5' adapter sequence
    """
    cdef cystr path, pathBS, pathCons
    from subprocess import check_call, CalledProcessError
    from itertools import chain
    from sys import stderr
    cfi = chain.from_iterable
    if(num_nucs < 0):
        raise UnsetRequiredParameter(
            "num_nucs required for DispatchParallelDMP.")
    elif(num_nucs == 0):
        raise ImproperArgumentError(
            "num_nucs shouldn't be set to 0. Otherwise, we'd just "
            "do it all in a single thread.")
    if(head < 0):
        raise UnsetRequiredParameter(
            "head required for DispatchParallelDMP.")
    if(threads < 0):
        raise UnsetRequiredParameter(
            "threads required for DispatchParallelDMP.")
    elif(threads == 1):
        raise ImproperArgumentError(
            "threads shouldn't be set to 0. Otherwise, we'd just "
            "do it all in a single thread and avoid this Popen rat's nest.")
    # Split the fastqs by the start of the barcode sequence
    fqPairList = PairedShadeSplitter(fq1, fq2, indexFq=indexFq,
                                     head=head, num_nucs=num_nucs)
    # Create the PopenDispatcher - this submits both foreground and
    # background jobs.
    Dispatcher = GetParallelDMPPopen(fqPairList, sortMem=sortMem,
                                     threads=threads, head=head,
                                     overlapLen=overlapLen,
                                     p3Seq=p3Seq, p5Seq=p5Seq,
                                     cap=cap, cap_quality=cap_quality)
    Dispatcher.daemon()
    outfqpaths = [(i.split(",")[0], i.split(",")[1]) for
                  i in Dispatcher.outstrs.itervalues()]
    outfq1s = [i[0] for i in outfqpaths]
    outfq2s = [i[1] for i in outfqpaths]
    del outfqpaths
    outfq1 = TrimExt(fq1) + ".dmp.merged.fastq"
    outfq2 = TrimExt(fq2) + ".dmp.merged.fastq"
    catStr1 = " ".join(["cat"] + outfq1s + [">", outfq1])
    pl("About to clean up my R1 temporary files with command %s." % catStr1)
    check_call(catStr1, shell=True, executable="/bin/bash")
    check_call(shlex.split("rm " + " ".join(outfq1s)))
    catStr2 = " ".join(["cat"] + outfq2s + [">", outfq2])
    pl("About to clean up my R2 temporary files with command %s." % catStr2)
    check_call(catStr2, shell=True, executable="/bin/bash")
    check_call(shlex.split("rm " + " ".join(outfq2s)))
    for path in cfi(fqPairList):
        pathBS = path.replace(".fastq", ".BS.fastq")
        pathCons = path.replace(".fastq", ".BS.cons.fastq")
        try:
            stderr.write("Now calling 'rm %s %s %s'\n" % path, pathBS, pathCons)
            check_call(["rm", path, pathBS, pathCons])
        except CalledProcessError:
            raise Exception("Path attempting to remove: %s\n" % path)
    return outfq1, outfq2

REMOVE_NS = maketrans("N", "A")

def PairedShadeSplitter(cystr fq1, cystr fq2, cystr indexFq="default",
                        int head=-1, int num_nucs=-1):
    """
    :param [cystr/arg] fq1 - path to read 1 fastq
    :param [cystr/arg] fq2 - path to read 2 fastq
    :param [cystr/kwarg/"default"] indexFq - path to index fastq
    :param [int/kwarg/2] head - number of bases each from reads 1
    and 2 with which to salt the barcodes.
    :param [object/nkwarg/-1] num_nucs - number of nucleotides at the
    beginning of a barcode to include in creating the output handles.
    """
    # Imports
    from utilBMF.HTSUtils import nci
    #  C declarations
    cdef pFastqProxy_t read1
    cdef pFastqProxy_t read2
    cdef pysam.cfaidx.FastqProxy indexRead
    cdef list bcKeys
    cdef dict BarcodeHandleDict1, BarcodeHandleDict2
    cdef int hpLimit
    if(head < 0):
        raise UnsetRequiredParameter(
            "PairedShadeSplitting requires that head be set.")
    if(num_nucs < 0):
        raise UnsetRequiredParameter("num_nucs must be set")
    elif(num_nucs > 2):
        raise ImproperArgumentError("num_nucs is limited to 2")
    numHandleSets = 4 ** num_nucs
    bcKeys = ["A" * (num_nucs - len(nci(i))) + nci(i) for
              i in range(numHandleSets)]
    if(indexFq is None):
        raise UnsetRequiredParameter(
            "indexFq required for PairedShadeSplitter.")
    base_outfq1 = TrimExt(fq1).replace(".fastq",
                                       "").split("/")[-1] + ".shaded."
    base_outfq2 = TrimExt(fq2).replace(".fastq",
                                       "").split("/")[-1] + ".shaded."
    BarcodeHandleDict1 = {key: open(base_outfq1 + key +
                                    ".fastq", "w") for key in bcKeys}
    BarcodeHandleDict2 = {key: open(base_outfq2 + key +
                                    ".fastq", "w") for key in bcKeys}
    inFq1 = pFastqFile(fq1)
    inFq2 = pFastqFile(fq2)
    ifn2 = inFq2.next
    inIndex = pysam.FastqFile(indexFq, persist=False)
    hpLimit = len(inIndex.next().sequence) * 3 // 4
    inIndex = pysam.FastqFile(indexFq, persist=False)
    ifin = inIndex.next
    numWritten = 0
    for read1 in inFq1:
        read2 = ifn2()
        indexRead = ifin()
        tempBar = (read1.sequence[1:head + 1] + indexRead.sequence +
                   read2.sequence[1:head + 1])
        bin = tempBar[:num_nucs]
        read1.comment = cMakeTagComment(tempBar, read1, hpLimit)
        read2.comment = cMakeTagComment(tempBar, read2, hpLimit)
        BarcodeHandleDict1[bin.translate(REMOVE_NS)].write(str(read1))
        BarcodeHandleDict2[bin.translate(REMOVE_NS)].write(str(read2))
    [outHandle.close() for outHandle in BarcodeHandleDict1.itervalues()]
    [outHandle.close() for outHandle in BarcodeHandleDict2.itervalues()]
    return zip([i.name for i in BarcodeHandleDict1.itervalues()],
               [j.name for j in BarcodeHandleDict2.itervalues()])


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
    inIndex = pysam.FastqFile(indexFq, persist=False)
    hpLimit = len(inIndex.next().sequence) * 3 // 4
    inIndex = pysam.FastqFile(indexFq, persist=False)
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
        tempBar = (read1.sequence[1:head + 1] + indexRead.sequence +
                   read2.sequence[1:head + 1])
        read1.comment = cMakeTagComment(tempBar, read1, hpLimit)
        read2.comment = cMakeTagComment(tempBar, read2, hpLimit)
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
        pIndexRead = pFastqProxy.fromFastqProxy(inIndex.next())
        pRead1 = pFastqProxy.fromFastqProxy(read1)
        if("N" in pRead1.pIndexRead.sequence):
            read1.comment += "|FP=0|BS=%" % (
                read1.sequence[1:head + 1] + pIndexRead.sequence)
            outFqHandle1.write(str(read1))
        else:
            read1.comment += "|FP=1|BS=%s" % (
                read1.sequence[1:head + 1] + pIndexRead.sequence)
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


cdef inline cystr CappedPyArr2QualStr(py_array qualArr, int8_t cap=38):
    """
    :param qualArr: [py_array/arg] Expects a 32-bit integer array
    :return: [cystr]
    """
    cdef size_t i
    cdef int8_t cap_plus_33 = cap + 33
    cdef py_array ret = array("B", [i + 33 if(i < cap) else cap_plus_33 for
                                    i in qualArr])
    return ret.tostring()


cdef inline cystr PyArr2QualStr(py_array qualArr):
    """
    :param qualArr: [py_array/arg] Expects a 32-bit integer array
    :return: [cystr]
    """
    cdef size_t i
    cdef py_array ret = array("B", [
        i + 33 if(i < 94) else 126 for i in qualArr])
    return ret.tostring()


cdef cystr cQualArr2QualStr(ndarray[int32_t, ndim=1] qualArr):
    """
    This is the "safe" way to convert ph2chr.
    """
    cdef int32_t tmpInt
    return array('B', [
        tmpInt if(tmpInt < 94) else
        93 for tmpInt in qualArr]).tostring().translate(PH2CHR_TRANS)


cdef inline cystr PyArr2FAString(py_array agreeArr):
    cdef int32_t i
    return "|FA=%s" % ",".join([str(i) for i in agreeArr])


cdef inline cystr PyArr2PVString(py_array qualArr):
    cdef int32_t i
    return "|PV=%s" % ",".join([str(i) for i in qualArr])


cdef cystr cQualArr2PVString(ndarray[int32_t, ndim=1] qualArr):
    return "|PV=%s" % ",".join(qualArr.astype(str))


cdef cystr cQualArr2FAString(ndarray[int32_t, ndim=1] qualArr):
    return "|FA=%s" % ",".join(qualArr.astype(str))


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


def pairedFastqConsolidate(fq1, fq2,
                           int SetSize=100, bint parallel=True,
                           bint cap_quality=False,
                           int8_t cap=38):
    """
    TODO: Unit test for this function.
    Also, it would be nice to do a groupby() that separates read 1 and read
    2 records so that it's more pythonic, but that's a hassle.
    """
    cdef cystr outFq1, outFq2, cStr
    cdef int checks
    pl("Now running pairedFastqConsolidate on {} and {}.".format(fq1, fq2))
    pl("(What that really means is that I'm running "
       "singleFastqConsolidate twice")
    if(not parallel):
        outFq1 = singleFastqConsolidate(fq1,
                                        SetSize=SetSize)
        outFq2 = singleFastqConsolidate(fq2,
                                        SetSize=SetSize)
    else:
        # Make background process command string
        cStr = ("python -c 'from MawCluster.BCFastq import singleFastqConsoli"
                "date;singleFastqConsolidate(\"%s\", SetSize=" % fq2 +
                "%s, cap_quality=%s, cap=%s);" % (SetSize, cap_quality, cap) +
                "import sys;sys.exit(0)'")
        # Submit
        PFC_Call = subprocess.Popen(cStr, stderr=None, shell=True,
                                    stdout=None, stdin=None, close_fds=True)
        # Run foreground job on other read
        outFq1 = singleFastqConsolidate(fq1, SetSize=SetSize,
                                        cap_quality=cap_quality, cap=cap)
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
    return outFq1, outFq2


def singleFqConsRescale(cystr fq,
                        cystr RescaleData,
                        int SetSize=100):
    """
    Takes the path to a python pickle for the rescaling array.
    :param fq:
    :param SetSize:
    :param RescaleData:
    :return:
    """
    cdef cystr outFq, bc4fq, ffq
    cdef pFastqFile_t inFq
    cdef list StringList, pFqPrxList
    cdef int numproc, TotalCount, MergedCount
    cdef py_array RescalingArray

    from sys import stderr

    pl("Now running singleFastqConsolidate on {}.".format(fq))
    RescalingArray = load(open(RescaleData, "rb"))
    outFq = TrimExt(fq) + ".cons.fastq"
    inFq = pFastqFile(fq)
    outputHandle = open(outFq, 'w')
    StringList = []
    numproc = 0
    TotalCount = 0
    ohw = outputHandle.write
    sla = StringList.append
    for bc4fq, fqRecGen in groupby(inFq, key=getBS):
        pFqPrxList = list(fqRecGen)
        ffq = cFastFisherFlattenRescale(pFqPrxList, bc4fq, RescalingArray)
        sla(ffq)
        numproc += 1
        TotalCount += len(pFqPrxList)
        MergedCount += 1
        if not (numproc % SetSize):
            ohw("".join(StringList))
            StringList = []
            sla = StringList.append
            continue
    ohw("".join(StringList))
    outputHandle.flush()
    inFq.close()
    outputHandle.close()
    if("SampleMetrics" in globals()):
        globals()['SampleMetrics']['TotalReadCount'] = TotalCount
        globals()['SampleMetrics']['MergedReadCount'] = MergedCount
    stderr.write("Consolidation a success for inFq: %s!\n" % fq)
    return outFq


def singleFastqConsolidate(cystr fq,
                           int SetSize=100,
                           bint cap_quality=False,
                           int8_t cap=38):
    cdef cystr outFq, bc4fq, ffq
    cdef pFastqFile_t inFq
    cdef list StringList, pFqPrxList
    cdef int numproc, TotalCount, MergedCount
    from sys import stderr
    outFq = TrimExt(fq) + ".cons.fastq"
    pl("Now running singleFastqConsolidate on {}.".format(fq))
    inFq = pFastqFile(fq)
    outputHandle = open(outFq, 'w')
    StringList = []
    numproc = 0
    TotalCount = 0
    ohw = outputHandle.write
    sla = StringList.append
    for bc4fq, fqRecGen in groupby(inFq, key=getBS):
        pFqPrxList = list(fqRecGen)
        ffq = cFastFisherFlattening(pFqPrxList, bc4fq,
                                    cap_quality=cap_quality, cap=cap)
        sla(ffq)
        numproc += 1
        TotalCount += len(pFqPrxList)
        MergedCount += 1
        if not (numproc % SetSize):
            ohw("".join(StringList))
            StringList = []
            sla = StringList.append
            continue
    ohw("".join(StringList))
    outputHandle.flush()
    inFq.close()
    outputHandle.close()
    if("SampleMetrics" in globals()):
        globals()['SampleMetrics']['TotalReadCount'] = TotalCount
        globals()['SampleMetrics']['MergedReadCount'] = MergedCount
    stderr.write("Consolidation a success for inFq: %s!\n" % fq)
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
            ew(str(pFastqProxy.fromFastqProxy(read)))
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
            ew(str(pFastqProxy.fromFastqProxy(read1)))
            ew(str(pFastqProxy.fromFastqProxy(read2)))
            continue
        if homing not in read2.sequence[bcLen:bcLen + HomingLen]:
            ew(str(pFastqProxy.fromFastqProxy(read1)))
            ew(str(pFastqProxy.fromFastqProxy(read2)))
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
    TrueFamDict = {x: x for x, y1 in histDict.iteritems() if y1 >= minFam}
    if(len(TrueFamDict) == 0):
        rescueHistDict = {}
    else:
        rescueHistDict = {h: y for y in TrueFamDict.iterkeys() for
                          h in hamming_cousins(y, n=n)}
    if(tmpFile is not None):
        check_call(["rm", tmpFile])
    return rescueHistDict, TrueFamDict


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
    cdef cystr tmpBS, saltedBS, tagStr, indexSeq
    cdef pFastqFile_t inHandle1, inHandle2, indexHandle
    cdef pFastqProxy_t rec1, rec2, index_read
    cdef int bLen, hpLimit
    if(outFq1 is None):
        outFq1 = TrimExt(inFq1, exclude="fastq") + ".rescued.shaded.fastq"
    if(outFq2 is None):
        outFq2 = TrimExt(inFq2, exclude="fastq") + ".rescued.shaded.fastq"
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
    hpLimit = bLen * 3 // 4  # Homopolymer limit
    indexHandle = pFastqFile(indexFq)
    ihn = indexHandle.next
    for rec1 in inHandle1:
        try:
            index_read = ihn()
            rec2 = ih2n()
        except StopIteration:
            raise Tim("Index fastq and read fastqs have different sizes. "
                      "Abort!")
        try:
            indexSeq = TrueFamDict[index_read.sequence]
            saltedBS = (rec1.sequence[1:head + 1] + indexSeq +
                        rec2.sequence[1:head + 1])
            rec1.comment = cMakeTagComment(saltedBS, rec1, hpLimit)
            rec2.comment = cMakeTagComment(saltedBS, rec2, hpLimit)
        except KeyError:
            pass
            try:
                indexSeq = rescueDict[index_read.sequence]
                saltedBS = (rec1.sequence[1:head + 1] +
                            indexSeq +
                            rec2.sequence[1:head + 1])
                rec1.comment = cMakeTagComment(saltedBS, rec1, hpLimit)
                rec2.comment = cMakeTagComment(saltedBS, rec2, hpLimit)
            except KeyError:
                # This isn't in a true family. Blech!
                saltedBS = (rec1.sequence[1:head + 1] + index_read.sequence +
                            rec2.sequence[1:head + 1])
                rec1.comment = cMakeTagComment(saltedBS, rec1, hpLimit)
                rec2.comment = cMakeTagComment(saltedBS, rec2, hpLimit)
        ohw1(str(rec1))
        ohw2(str(rec2))
    outHandle1.close()
    outHandle2.close()
    return outFq1, outFq2


cpdef cystr MakeTagComment(cystr saltedBS, pFastqProxy_t rec, int hpLimit):
    """
    Python-visible MakeTagComment.
    """
    return cMakeTagComment(saltedBS, rec, hpLimit)


cdef inline bint BarcodePasses(cystr barcode, int hpLimit):
    return not ("N" in barcode or "A" * hpLimit in barcode or
                "C" * hpLimit in barcode or
                "G" * hpLimit in barcode or "T" * hpLimit in barcode)
