#!/usr/bin/env python
from __future__ import division

import logging
import os.path
import subprocess
import shlex
import operator
from operator import attrgetter as oag
from operator import methodcaller as mc
from subprocess import check_output, check_call
import sys

import cython
cimport cython
import pysam
cimport pysam.cfaidx
cimport pysam.calignmentfile
import numpy as np
cimport numpy as np
from numpy import array as nparray
from numpy import append as npappend
from numpy import mean as nmean
from cytoolz import map as cmap

from utilBMF.HTSUtils import (ParseBed, printlog as pl, CoorSortAndIndexBam,
                              ThisIsMadness, PipedShellCall)
from MawCluster.PileupUtils import pPileupColumn

ctypedef np.longdouble_t dtype128_t
ctypedef np.int64_t dtypei64_t


"""
Contains functions and miscellania for QC metrics.
"""

@cython.locals(buffer=cython.long)
def ExtendBed(bedfile, buffer=100, outbed="default"):
    """
    """
    if outbed == "default":
        outbed = ".".join(bedfile.split(".")[0:-1] + ['extended',
                                                      'bed']).split("/")[-1]
    commandStr = ("cat %s | awk 'FS=OFS=\"\t\" {{print $1, $2 - " % bedfile +
                  "{0}, $3 + {0}, $4, $5, $6, $7, $8, $9, ".format(buffer) +
                  "$10, $11, $12}} > %s'" % outbed)
    PipedShellCall(outbed)
    return outbed


def coverageBed(inBAM, bed="default", outbed="default"):
    """
    Uses samtool bedcov tool.
    Requires samtools >= 1.0
    If Popen is true, instead of calling and waiting for the
    call to finish successfully, it opens a Popen instance which can be
    queried and returns a tuple of the Popen and the outbed. Use if you
    want this to go on in the background.
    """
    if(bed == "default"):
        raise ThisIsMadness("bed file required for calling "
                            "utilBMF.QC.coverageBed")
    if(outbed == "default"):
        outbed = ".".join(bed.split(".")[0:-1] + ["cov", "bed"])
    pl("Calculating coverage bed. bed: %s. inBAM: %s. " (bed, inBAM) +
       "outbed: %s" % outbed)
    commandStr = ("samtools bedcov %s %s | awk -t '{print " % (bed, inBAM) +
                  "$1, $2, $3, $4, $5, $6, $7 / ($3 - $2)} > %s" % ( outbed))
    PipedShellCall(commandStr)
    return outbed


@cython.locals(buffer=cython.long)
def FracOnTarget(inBAM, bed="default", buffer=100):
    """
    Calculates the fraction of mapped reads aligned to a target region.
    Do to "bleeding out" of regions of capture, there is an optional
    buffer option which counts reads just outside of the regions targeted
    as being "on target."
    """
    cdef cython.long numReadsTotal
    cdef cython.long numReadsOnTarget
    cdef cython.float fracOnTarget
    covBed = coverageBed(inBAM, bed=ExtendBed(bed, buffer=buffer))
    numReadsTotal = int(check_output(["samtools", "view", "-L", covBed,
                                      "-c", inBAM]).strip())
    numReadsOnTarget = int(check_output(
        "awk {print $NF} %s | paste -sd+ | bc" % covBed,
        shell=True).strip())
    fracOnTarget = numReadsOnTarget / (1. * numReadsTotal)
    return covBed, fracOnTarget


@cython.returns(cython.long)
def CountNumReads(inBAM):
    """
    Simply counts the number of reads in a BAM file.
    """
    return int(check_output(["samtools", "view", "-c", inBAM]).strip())


def InsertSizeArray(inBAM):
    """
    Returns an array of insert sizes for the sample.
    """
    cdef pysam.calignmentfile.AlignedSegment i
    inHandle = pysam.AlignmentFile(inBAM, "rb")
    return [abs(i.tlen) for i in inHandle]


@cython.locals(min=cython.long, onTargetBuffer=cython.long)
def GetAllQCMetrics(inBAM, bedfile="default", onTargetBuffer=100,
                    minFM=2):
    cdef cython.long TotalReads
    cdef cython.long MappedReads
    cdef cython.long UnmappedReads
    cdef cython.float fracOnTarget
    cdef cython.float meanInsert
    cdef cython.long MappedFamReads
    cdef cython.long MappedSingletonReads
    cdef np.ndarray[dtypei64_t, ndim = 1] allInserts
    outfile = ".".join(inBAM.split(".")[0:-1] + ["qc", "txt"])
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
    outHandle.write("#" +": ".join(["BAM", inBAM]) + "\n")
    outHandle.write("#" +": ".join(["bed", bedfile]) + "\n")
    outHandle.write("#" +": ".join(["Expanded bed", extendedBed]) + "\n")
    outHandle.write("#" +": ".join(["Target Coverage Buffer",
                                    str(onTargetBuffer)]) + "\n")
    outHandle.write("#" +": ".join(["Minimum Family Coverage", str(minFM)]) + "\n")
    popenInstance = coverageBed(inBAM, bed=bedfile, )
    outHandle.write("#" + ": ".join(["Output coverage bed", bedfile]))
    TotalReads = CountNumReads(inBAM)
    resultsDict["TotalReads"] = TotalReads
    MappedReads = 0
    ReadsOnTarget = 0
    ReadsOffTarget = 0
    MappedFamReads = 0
    if(os.path.isfile(inBAM + ".bai") is False):
        check_call(["samtools", "index", inBAM])
        if(os.path.isfile(inBAM + ".bai") is False):
            pl("Input BAM was not coordinate sorted - doing so now.")
            inBAM = CoorSortAndIndexBam(inBAM)
    covBed, fracOnTarget = coverageBed(inBAM, bed=extendedBed)
    resultsDict["covbed"] = covBed
    resultsDict["fracOnTarget"] = fracOnTarget
    inHandle = pysam.AlignmentFile(inBAM, "rb")
    pileupIterator = inHandle.pileup(max_depth=200000, multiple_iterators=True)
    puIt = pileupIterator.next
    allInserts = nparray([], dtype=np.int64)
    allFMs = nparray([], dtype=np.int64)
    while True:
        try:
            p = pPileupColumn(puIt())
            pileups = p.pileups
            MappedReads += p.n
            allInserts = npappend(allInserts, list(cmap(
                oag("alignment.template_length", pileups))))
            allFMs = npappend(allFMs, list(cmap(
                mc("opt", "FM"), list(cmap(oag("alignment"), pileups)))))
            MappedFamReads += len([i for i in pileups if
                                   i.alignment.opt("FM") >= minFM])
        except StopIteration:
            pl("Finished iterating through the BAM file. Exiting loop.")
            break
    UnmappedReads = TotalReads - MappedReads
    resultsDict["MappedFamReads"] = MappedFamReads
    resultsDict["UnmappedReads"] = UnmappedReads
    MappedSingletonReads = MappedReads - MappedFamReads
    allInserts = allInserts[allInserts > 0]  # Remove tlen's of 0 from mean.
    meanInsert = nmean(allInserts)
    resultsDict["MeanInsert"] = meanInsert
    resultsDict["MappedSingletonReads"] = MappedSingletonReads
    outHandle.write("Mean Insert Size=%s\n" % meanInsert + "\n")
    outHandle.write("Unmapped Reads=%s\n" % UnmappedReads + "\n")
    outHandle.write("Mapped Reads=%s\n" % MappedReads + "\n")
    outHandle.write("Fraction On Target=%s\n" % fracOnTarget)
    outHandle.write(
        "Mapped Reads With Family Size g/e %s=%s\n" % (minFM, MappedReads))
    outHandle.write(
        "Mapped Reads With Family Size less than "
        "%s=%s\n" % (minFM, MappedSingletonReads))
    outHandle.close()
    return resultsDict


def GetFamSizeStats(inFq, outfile=sys.stdout):
    """
    Calculates family size and library diversity from a consolidated
    fastq file.
    """
    if(isinstance(outfile, str)):
        outfile = open(outfile, "w")
    cdef pysam.cfaidx.FastqProxy read
    cdef cython.long numFam = 0
    cdef cython.long numSing = 0
    cdef cython.long sumFam = 0
    cdef cython.long sumAll = 0
    cdef cython.long famS
    cdef cython.float MeanFamAll
    cdef cython.float MeanRealFam
    cdef pysam.cfaidx.FastqFile FqHandle
    FqHandle = pysam.FastqFile(inFq)
    for read in FqHandle:
        famS = int(read.comment.split(" ")[3].split("=")[1])
        if(famS > 1):
            numFam += 1
            sumFam += famS
        else:
            numSing += 1
        sumAll += famS
    MeanFamAll = sumAll / (1. * (numFam + numSing))
    MeanRealFam = sumFam / (1. * numFam)
    string = ("numFam:%s. NumSing:%s. " % (numFam, numSing) +
              "MeanFamAll:%s. MeanRealFam: %s" % (MeanFamAll, MeanRealFam))
    outfile.write(string + "\n")
    outfile.close()
    return numFam, numSing, MeanFamAll, MeanRealFam
