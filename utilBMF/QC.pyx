#!/usr/bin/env python
from __future__ import division

import logging
import os.path
import subprocess
import shlex
import operator
from operator import attrgetter as oag
from operator import methodcaller as mc

import cython
cimport cython
import pysam
cimport pysam.cfaidx
import numpy as np
cimport numpy as np
from numpy import array as nparray
from numpy import append as npappend
from numpy import mean as nmean

from utilBMF.HTSUtils import (ParseBed, printlog as pl, CoorSortAndIndexBam,
                              ThisIsMadness, PipedShellCall)
from MawCluster.PileupUtils import pPileupColumn

ctypedef np.longdouble_t dtype128_t
ctypedef np.int64_t dtypei64_t


"""
Contains functions and miscellania for QC metrics.
"""


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
    


@cython.locals(Popen=cython.bint)
def coverageBed(inBAM, bed="default", outbed="default",
                Popen=False):
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
    commandStr = "samtools bedcov %s %s > %s" % (bed, inBAM)
    if(Popen is False):
        PipedShellCall(commandStr)
        return outbed
    else:
        outHandle = open(outbed, "w")
        return subprocess.Popen(shlex.split(
            "samtools bedcov %s %s" % (bed, inBAM)),
                                stdout=outHandle), outbed


@cython.locals(min=cython.long)
def meanFamSize(inBAM, min=2):
    """
    Calculates the mean family size for reads with a family size greater
    than or equal to the minimum. Default: 2
    """
    pass


@cython.locals(buffer=cython.long)
def FracOnTarget(inBAM, bed="default", buffer=100):
    """
    Calculates the fraction of mapped reads aligned to a target region.
    Do to "bleeding out" of regions of capture, there is an optional
    buffer option which counts reads just outside of the regions targeted
    as being "on target."
    """
    pass


@cython.returns(cython.long)
def NumUniqueReads(inBAM):
    """
    returns the number of reads with family size of 1
    """
    pass


def MeanInsertSize(inBAM):
    """
    Returns the mean insert size for the sample.
    """
    pass



@cython.locals(min=cython.long, onTargetBuffer=cython.long)
def GetAllQCMetrics(inBAM, bedfile="default", onTargetBuffer=100,
                    min=2):
    cdef cython.long TotalReads
    cdef cython.long MergedReads
    cdef cython.long MappedReads
    cdef cython.long ReadsOnTarget
    cdef cython.long ReadsOffTarget
    cdef np.ndarray[dtypei64_t, ndim = 1] allInserts
    outfile = ".".join(inBAM.split(".")[0:-1] + ["qc", "txt"])
    outHandle = open(outfile, "w")
    pl("GetAllQCMetrics running")
    if(bedfile == "default"):
        raise ThisIsMadness("bedfile must be set of QC metrics!")
    if(onTargetBuffer > 0):
        extendedBed = ExtendBed(bedfile)
    else:
        extendedBed = bedfile
    bed = ParseBed(extendedBed)
    outHandle.write("##Parameters\n")
    outHandle.write("#" +": ".join(["BAM", inBAM]) + "\n")
    outHandle.write("#" +": ".join(["bed", bedfile]) + "\n")
    outHandle.write("#" +": ".join(["Expanded bed", extendedBed]) + "\n")
    outHandle.write("#" +": ".join(["Bed expansion",
                                    str(onTargetBuffer)]) + "\n")
    outHandle.write("#" +": ".join(["Target Coverage Buffer",
                                    str(onTargetBuffer)]) + "\n")
    outHandle.write("#" +": ".join(["Minimum Family Coverage", str(min)]) + "\n")
    popenInstance = coverageBed(inBAM, bed=bedfile, )
    outHandle.write("#" + ": ".join(["Output coverage bed", ]))
    TotalReads = 0
    MappedReads = 0
    ReadsOnTarget = 0
    ReadsOffTarget = 0
    if(os.path.isfile(inBAM + ".bai") is False):
        subprocess.check_call(["samtools", "index", inBAM])
        if(os.path.isfile(inBAM + ".bai") is False):
            inBAM = CoorSortAndIndexBam(inBAM)
    coverageBedHandle, covBed = coverageBed(inBAM,
                                            bed=extendedBed, Popen=True)
    inHandle = pysam.AlignmentFile(inBAM, "rb")
    pileupIterator = inHandle.pileup(max_depth=200000, multiple_iterators=True)
    firstIter = True
    while True:
        try:
            p = pPileupColumn(next(pileupIterator))
            pileups = p.pileups
            if(firstIter is True):
                allInserts = nparray(map(oag(
                    "template_length", map(
                        oag("alignment"), pileups))),
                               dtype=np.int64)
            else:
                allInserts = npappend(allInserts, map(
                    oag("template_length", map(
                        oag("alignment"), pileups))))
            if(firstIter is True):
                allFMs = nparray(map(mc("opt", "FM"), map(
                        oag("alignment"), pileups)), dtype=np.int64)
            else:
                allFMs = npappend(allFMs, map(
                    mc("opt", "FM"), map(oag("alignment"), pileups)))
            firstIter = False
            
        except StopIteration:
            pl("Finished iterations. (CalcWithoutBedCoverage)")
    allInserts = allInserts[allInserts > 0]  # Remove tlen's of 0 from mean.
    outHandle.write("mean insert size=%s\n" % (nmean(allInserts)))


def GetFamSizeStats(inFq):
    """
    Calculates family size and library diversity from a consolidated
    fastq file.
    """
    cdef pysam.cfaidx.FastqProxy read
    cdef cython.long numFam
    cdef cython.long numSing
    cdef cython.long sumFam
    cdef cython.long sumAll
    cdef cython.long famS
    cdef cython.float MeanFamAll
    cdef cython.float MeanRealFam
    cdef pysam.cfaidx.FastqFile FqHandle
    FqHandle = pysam.FastqFile(inFq)
    numFam = 0
    numSing = 0
    sumFam = 0
    sumAll = 0
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
    print("numFam:%s. NumSing:%s. MeanFamAll:%s. MeanRealFam: %s" % (
        numFam, numSing,
        MeanFamAll,
        MeanRealFam))
    return {name: eval(name) for name in ["numFam", "numSing", "MeanFamAll",
                                          "MeanRealFam"]}
