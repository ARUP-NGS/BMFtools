# cython: c_string_type=str, c_string_encoding=ascii
from __future__ import division
from operator import attrgetter as oag
from operator import methodcaller as mc
from subprocess import check_output, check_call
from numpy import (append as npappend,
                   mean as nmean, max as nmax)
from cytoolz import map as cmap
from .HTSUtils import (ParseBed, printlog as pl, CoorSortAndIndexBam,
                       pFastqFile, cyStdFlt, cyStdInt, TrimExt)
from MawCluster.BCFastq import GetDescriptionTagDict as descDict
from .ErrorHandling import ThisIsMadness, MissingExternalTool
from MawCluster.PileupUtils import pPileupColumn
from os import path as ospath
import cython
import logging
import numpy as np
import operator
import pysam
import shlex
import subprocess
import sys


"""
Contains functions and miscellania for QC metrics.
"""


@cython.locals(buffer=int)
def ExtendBed(bedfile, buffer=100, outbed="default"):
    """
    """
    if outbed == "default":
        outbed = ".".join(bedfile.split(".")[0:-1] + ['extended',
                                                      'bed']).split("/")[-1]
    commandStr = ("cat %s | awk 'FS=OFS=\"\t\" {{print $1, $2 - " % bedfile +
                  "{0}, $3 + {0}, $4, $5, $6".format(buffer) +
                  "}}' > %s" % outbed)
    pl("ExtendBed command string: %s" % commandStr)
    check_call(commandStr, shell=True)
    return outbed


def FastDOCBed(inBAM, bed="default", outbed="default",
               cystr FastDOCPath="default", int threads=8,
               int readLen=-1, int minMQ=-1):
    """
    Uses samtool bedcov tool.
    Requires samtools >= 1.
    If Popen is true, instead of calling and waiting for the
    call to finish successfully, it opens a Popen instance which can be
    queried and returns a tuple of the Popen and the outbed. Use if you
    want this to go on in the background.
    """
    cdef double fracOnTarget, numReadsOnTarget, numBasesOnTarget
    cdef int numReadsTotal
    if(readLen < 0):
        pl("readLen not set, inferring from inBAM %s." % inBAM)
        readLen = len(pysam.AlignmentFile(inBAM, "rb").next().seq)
        pl("readLen inferred to be %s" % readLen)
    if(bed == "default"):
        raise ThisIsMadness("bed file required for calling "
                            "utilBMF.QC.FastDOCBed")
    if(outbed == "default"):
        outbed = TrimExt(bed) + "cov.bed"
    if(minMQ < 0):
        pl("minMQ was not set for FastDOCBed. Defaulting to 1.")
        minMQ = 1
    if(FastDOCPath == "default"):
        try:
            print("Getting Chapman object")
            Chapman = globals()['Chapman']
            print("Getting FastDOCPath")
            FastDOCPath = Chapman['FastDOCPath']
        except KeyError:
            print("globals: %s" % globals())
            raise MissingExternalTool("FastDOCPath must be set to call FastDOC!")
    commandStr = ("java -jar %s %s %s " % (FastDOCPath, bed, inBAM) +
                  " --threads %s --format bed -m %s > " % (threads, minMQ) +
                  outbed)
    pl("Calculating coverage bed. bed: %s. inBAM: %s. " % (bed, inBAM) +
       "outbed: %s. Command string: %s" % (outbed, commandStr))
    check_call(commandStr, shell=True)
    commandStr = ("awk 'FS=OFS=\"\t\" {{print ($3 - $2) * $4}}' "
                  "%s | grep -v '^$'" % outbed)
    numBasesOnTarget = sum(float(i) for i in check_output(
        commandStr.replace("\t", "\\t"), shell=True).split("\n") if
                           i != "")
    pl("Number of bases on target: %s" % numBasesOnTarget)
    numReadsOnTarget = numBasesOnTarget / readLen
    numReadsTotal = CountNumReads(inBAM)
    fracOnTarget = numReadsOnTarget / numReadsTotal
    sys.stderr.write("Fraction on target: %s\n" % fracOnTarget)
    return outbed, fracOnTarget


@cython.locals(buffer=int)
def FracOnTarget(inBAM, bed="default", buffer=20,
                 cystr FastDOCPath="default", int threads=8):
    """
    Calculates the fraction of mapped reads aligned to a target region.
    Do to "bleeding out" of regions of capture, there is an optional
    buffer option which counts reads just outside of the regions targeted
    as being "on target."
    """
    cdef int numReadsTotal
    cdef int numReadsOnTarget
    cdef double fracOnTarget
    covBed, fracOnTarget = FastDOCBed(inBAM,
                                      bed=ExtendBed(bed, buffer=buffer),
                                      FastDOCPath=FastDOCPath,
                                      threads=threads)
    return FracOnTarget


@cython.returns(int)
def CountNumReads(inBAM):
    """
    Simply counts the number of reads in a BAM file.
    """
    return int(check_output(["samtools", "view", "-c", inBAM]).strip())


cdef ndarray[np.float64_t, ndim=1] InsertSizeArray_(
        pysam.calignmentfile.AlignmentFile handle):
    cdef pysam.calignmentfile.AlignedSegment i
    return np.absolute(np.array([i.tlen for i in handle], dtype=np.float64))


@cython.returns(ndarray)
def InsertSizeArray(cystr inBAM):
    """
    Returns an array of insert sizes for the sample.
    """
    cdef pysam.calignmentfile.AlignedSegment i
    cdef pysam.calignmentfile.AlignmentFile inHandle
    inHandle = pysam.AlignmentFile(inBAM, "rb")
    return InsertSizeArray_(inHandle)


def GetAllQCMetrics(inBAM, bedfile="default", int onTargetBuffer=20,
                    int minFM=2, FastDOCPath="default", int minMQ=-1):
    cdef int TotalReads, FM, MappedReads, UnmappedReads
    cdef double fracOnTarget, stdInsert, meanInsert
    cdef int MappedFamReads, maxInsert
    cdef int MappedSingletonReads
    outfile = TrimExt(inBAM) + ".qc.txt"
    outHandle = open(outfile, "w")
    pl("GetAllQCMetrics running")
    if(bedfile == "default"):
        raise ThisIsMadness("bedfile must be set of QC metrics!")
    if(onTargetBuffer > 0):
        extendedBed = ExtendBed(bedfile, buffer=onTargetBuffer)
    else:
        extendedBed = bedfile
    bed = ParseBed(extendedBed)
    resultsDict = {}
    outHandle.write("##Parameters\n")
    outHandle.write("#" + ": ".join(["BAM", inBAM]) + "\n")
    outHandle.write("#" + ": ".join(["bed", bedfile]) + "\n")
    outHandle.write("#" + ": ".join(["Expanded bed", extendedBed]) + "\n")
    outHandle.write("#" + ": ".join(["Target Coverage Buffer",
                                    str(onTargetBuffer)]) + "\n")
    outHandle.write("#" + ": ".join(["Minimum Family Coverage",
                                     str(minFM)]) + "\n")
    outHandle.write("#" + ": ".join(["Output coverage bed", bedfile]))
    TotalReads = CountNumReads(inBAM)
    resultsDict["TotalReads"] = TotalReads
    MappedReads = 0
    ReadsOnTarget = 0
    ReadsOffTarget = 0
    MappedFamReads = 0
    covBed, fracOnTarget = FastDOCBed(inBAM, bed=extendedBed,
                                      FastDOCPath=FastDOCPath,
                                      minMQ=minMQ)
    if(ospath.isfile(inBAM + ".bai") is False):
        check_call(["samtools", "index", inBAM])
        sys.stderr.write("Fraction on target: %s" % fracOnTarget)
        if(ospath.isfile(inBAM + ".bai") is False):
            pl("Input BAM was not coordinate sorted - doing so now.")
            inBAM = CoorSortAndIndexBam(inBAM)
    resultsDict["covbed"] = covBed
    resultsDict["fracOnTarget"] = fracOnTarget
    print("About to iterate through alignment file to get more QC metrics.")
    inHandle = pysam.AlignmentFile(inBAM, "rb")
    ihIterator = inHandle.next
    allInserts = []
    while True:
        try:
            p = ihIterator()
            MappedReads += 1
            allInserts.append(p.template_length)
            FM = p.opt("FM")
            if(FM >= 2):
                MappedFamReads += 1
        except StopIteration:
            pl("Finished iterating through the BAM file. Exiting loop.")
            break
    allInserts = np.array(allInserts, dtype=np.int64)
    UnmappedReads = TotalReads - MappedReads
    resultsDict["MappedFamReads"] = MappedFamReads
    resultsDict["UnmappedReads"] = UnmappedReads
    MappedSingletonReads = MappedReads - MappedFamReads
    allInserts = allInserts[allInserts > 0]  # Remove tlen's of 0 from mean.
    meanInsert = nmean(allInserts[allInserts < 5000])
    stdInsert = cyStdInt(allInserts[allInserts < 5000])
    maxInsert = nmax(allInserts)
    # Don't include huge inserts in mean calculation.
    resultsDict["MeanInsert"] = meanInsert
    resultsDict["MaxInsert"] = maxInsert
    resultsDict["MappedSingletonReads"] = MappedSingletonReads
    outHandle.write("Mean Insert Size=%s\n" % meanInsert)
    outHandle.write("Max Insert Size=%s\n" % maxInsert)
    outHandle.write("Std Dev of Insert Size (< 5000): %s" % stdInsert)
    outHandle.write("Unmapped Reads=%s\n" % UnmappedReads)
    outHandle.write("Mapped Reads=%s\n" % MappedReads)
    outHandle.write("Fraction On Target=%s\n" % fracOnTarget)
    outHandle.write(
        "Mapped Reads With Family Size g/e %s=%s\n" % (minFM, MappedReads))
    outHandle.write(
        "Mapped Reads With Family Size less than "
        "%s=%s\n" % (minFM, MappedSingletonReads))
    outHandle.close()
    return resultsDict


cdef ndarray[double] GetFamSizeStats_(pFastqFile_t FqHandle):
    """
    Calculates family size and library diversity from a consolidated
    fastq file.
    """
    cdef int famS, numFam, numSing, sumFam, sumAll
    cdef pFastqProxy_t read
    cdef double MeanFamAll, MeanRealFam
    cdef cystr key, value
    cdef dict rDict
    numFam = 0
    numSing = 0
    sumFam = 0
    sumAll = 0
    for read in FqHandle:
        rDict = descDict(read.comment)
        famS = int(rDict["FM"])
        if(famS > 1):
            numFam += 1
            sumFam += famS
        else:
            numSing += 1
        sumAll += famS
    MeanFamAll = sumAll / (1. * (numFam + numSing))
    MeanRealFam = sumFam / (1. * numFam)
    return np.array([float(numFam), float(numSing), MeanFamAll,
                     MeanRealFam], dtype=np.double)


def GetFamSizeStats(cystr inFq, outfile=sys.stdout):
    """
    Calculates family size and library diversity from a consolidated
    fastq file.
    """
    cdef int numFam, numSing
    cdef pysam.cfaidx.FastqProxy read
    cdef double MeanFamAll, MeanRealFam
    cdef pFastqFile_t FqHandle
    cdef ndarray[double, ndim=1] results
    if(outfile == sys.stdout):
        outStr = "stdout"
    elif(isinstance(outfile, str)):
        outStr = outfile
        outfile = open(outfile, "w")
    else:
        outStr = repr(outfile)
    pl("Getting family stats for %s, outputting to %s" % (inFq, outStr))
    FqHandle = pFastqFile(inFq)
    results = GetFamSizeStats_(FqHandle)
    numFam = int(results[0])
    numSing = int(results[1])
    MeanFamAll = results[2]
    MeanRealFam = results[3]
    string = ("numFam:%s\nNumSing:%s\n" % (numFam, numSing) +
              "MeanFamAll:%s\nMeanRealFam: %s" % (MeanFamAll, MeanRealFam))
    outfile.write(string + "\n")
    outfile.close()
    return numFam, numSing, MeanFamAll, MeanRealFam
