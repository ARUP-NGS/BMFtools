# cython: c_string_type=str, c_string_encoding=ascii
# cython: cdivision=True

from __future__ import division
import subprocess
from os import path as ospath
import shlex
import logging
import operator
import warnings
from operator import attrgetter as oag
from operator import itemgetter as oig
from operator import div as odiv
from itertools import ifilterfalse as iff
from array import array
try:
    import re2 as re
except ImportError:
    import re

import numpy as np
from numpy import (mean as nmean, max as nmax, sum as nsum,
                   std as nstd)
import pysam
from pysam.calignmentfile import PileupRead as cPileupRead
import cython
from cytoolz import (frequencies as cyfreq,
                     partition as ctpartition)
from itertools import groupby
from utilBMF.ErrorHandling import ThisIsMadness, AbortMission
from utilBMF.HTSUtils import (ReadPair, printlog as pl, pPileupRead,
                              PileupReadPair, ssStringFromRead)
from utilBMF import HTSUtils
import utilBMF

#  Pre-defining some itemgetters and attrgetters
oig1 = oig(1)
oig0 = oig(0)
oagbc = oag("BaseCall")
oagir = oag("is_reverse")
oagread = oag("read")
oagqp = oag("query_position")
oagmq = oag("MQ")
oagbq = oag("BQ")
oagname = oag("name")
oagdisc = oag("discordant")
oagqn = oag("query_name")
oagal = oag("alignment")
oagname = oag("name")


@cython.returns(list)
def GetDiscReadNames(pPileupColumn_t pPC, oagname=oagname):
    cdef list returnList, gList
    cdef cystr name
    returnList = []
    for name, group in groupby(sorted(pPC.pileups, key=oagname), key=oagname):
        gList = list(group)
        if(len(gList) == 2 and gList[0].BaseCall != gList[1].BaseCall):
            returnList.append(name)
    return returnList


@cython.locals(paired=bint)
@cython.returns(list)
def GetDiscordantReadPairs(pPileupColumn_t pPileupColObj,
                           object oagname=oagname,
                           object oagdisc=oagdisc,
                           object ctpartition=ctpartition,
                           object oagal=oagal, object oagqn=oagqn):
    """
    Takes a pPileupColumn object (python PileupColumn) as input
    and returns a list of PileupReadPair objects.
    """
    cdef dict ReadNameCounter
    cdef list reads, readnames, readpairs, pileups
    cdef pPileupRead_t read
    cdef tuple i
    pileups = pPileupColObj.pileups
    ReadNameCounter = cyfreq(map(oagqn, map(oagal, pileups)))
    readnames = [i[0] for i in ReadNameCounter.iteritems() if i[1] == 2]
    reads = sorted([read for read in pileups if read.name in readnames],
                   key=oagname)
    readpairs = map(PileupReadPair, ctpartition(2, reads))
    return [pair for pair in readpairs if pair.discordant]


cdef class pPileupColumn:
    """
    Python container for the PileupColumn proxy in pysam.
    """
    def __cinit__(self, pysam.calignmentfile.PileupColumn PileupColumn):
        cdef PileupRead_t p
        self.pileups = [pPileupRead(p) for p in PileupColumn.pileups if not
                        p.is_del and not p.is_refskip]
        self.nsegments = len(self.pileups)
        self.reference_id = PileupColumn.reference_id
        self.reference_pos = PileupColumn.reference_pos


"""
Contains various utilities for working with pileup generators
and preparing for variant calls.
"""


cdef class PRInfo:

    """
    Created from a pysam.PileupRead object or its python equivalent,
    a pPileupRead.
    Holds family size, SV tags, base quality, mapping quality, and base.
    If any of the "finicky" fields aren't filled (e.g., if BAMs are not
    produced using BMFTools), they are set to None.
    Check to see if an attribute is None first if you'd like
    to save yourself some headaches.
    """

    cpdef object opt(self, cystr arg):
        return self.read.opt(arg)

    def __init__(self, pPileupRead_t PileupRead):
        cdef pysam.calignmentfile.AlignedSegment alignment
        alignment = PileupRead.alignment
        aopt = alignment.opt
        self.Pass = True
        self.read = alignment
        self.AF = aopt("AF")
        self.ND = aopt("ND")
        self.NF = aopt("NF")
        self.FM = aopt("FM")
        self.BQ = PileupRead.BQ
        self.FA = PileupRead.FA
        self.MQ = alignment.mapq
        try:
            self.SVPass = TestSVTags(aopt("SV"))
        except KeyError:
            self.SVPass = True
            warnings.warn(
                "Note: SV tag not present. All reads "
                "will pass the SV filter.", RuntimeWarning)
        if(alignment.is_qcfail or aopt("FP") == 0):
            self.Pass = False
        self.query_name = PileupRead.name
        self.BaseCall = alignment.query_sequence[
            PileupRead.query_position]
        if(self.BaseCall == "N"):
            self.Pass = False
        self.is_reverse = alignment.is_reverse
        self.is_proper_pair = alignment.is_proper_pair
        self.ssString = ssStringFromRead(alignment)
        self.query_position = PileupRead.query_position
        self.FractionAgreed = self.FA / (1. * self.FM)
        if(PileupRead.MBQ > 0):
            self.PVFrac = self.BQ * 1. / PileupRead.MBQ
        else:
            self.PVFrac = 0.


@cython.returns(cystr)
def is_reverse_to_str(bint boolean):
    if(boolean):
        return "reverse"
    if(boolean is False):
        return "forward"
    return "unmapped"


cdef class AlleleAggregateInfo:

    """
    Class which holds summary information for a given alt allele
    at a specific position. Meant to be used as part of the PCInfo
    class as values in the VariantDict.
    All alt alleles in this set of reads should be identical.
    recList must be a list of PRInfo objects.

    """
    def __init__(self, list recList, cystr consensus="default",
                 int mergedSize=-1, int totalSize=-1,
                 int minMQ=0,
                 int minBQ=0,
                 cystr contig="default",
                 int pos=-1,
                 int DOC=-1,
                 int DOCTotal=-1,
                 int NUMALT=-1,
                 float AABPSD=-1., float AAMBP=-1.,
                 float minFracAgreed=0.0, int minFA=0,
                 float minPVFrac=0.0, int FSR=-1,
                 object oagir=oagir, object oagqp=oagqp, object oagbq=oagbq,
                 object oagmq=oagmq, object nparray=np.array):
        cdef ndarray NFList
        cdef ndarray[np.float64_t, ndim=1] query_positions, PVFArray
        cdef int lenR
        cdef PRInfo_t rec
        cdef tuple i
        if(consensus == "default"):
            raise ThisIsMadness("A consensus nucleotide must be provided.")
        if(NUMALT < 0):
            raise ThisIsMadness("Number of alternate alleles required.")
        if(DOC < 0):
            raise ThisIsMadness("Full depth of coverage must be provided.")
        if(contig == "default"):
            raise ThisIsMadness("A contig must be provided (string).")
        if(pos < 0):
            raise ThisIsMadness("A position (0-based) must be provided (int).")
        if(mergedSize < 0):
            raise ThisIsMadness(("mergedSize must be provided: number of "
                                 "PRInfo at given position."))
        if(totalSize < 0):
            raise ThisIsMadness(("totalSize must be provided: number of "
                                 "PRInfo (pre-dmp) at given position."))
        self.NUMALT = NUMALT
        self.DOC = DOC
        self.DOCTotal = DOCTotal
        self.contig = contig
        self.pos = pos
        self.recList = [rec for rec in recList if
                        rec.MQ >= minMQ and rec.BQ >= minBQ and
                        rec.FractionAgreed >= minFracAgreed and
                        rec.FA >= minFA]
        # Check that all alt alleles are identical
        lenR = len(self.recList)
        self.len = lenR
        # Total Number of Differences
        if(lenR != 0):
            self.ALT = self.recList[0].BaseCall
            self.TND = sum(map(oag("ND"), self.recList))
            NFList = np.array(map(oag("NF"), self.recList))
            self.MNF = nmean(NFList)
            self.maxNF = nmax(NFList)
            self.NFSD = nstd(NFList)
        else:
            self.ALT = "N"
            self.TND = -1
            NFList = np.array([])
            self.MNF = -1.
            self.maxNF = -1.
            self.NFSD = -1.
        self.TotalReads = sum([rec.FM for rec in self.recList])
        self.MergedReads = lenR
        self.ReverseMergedReads = sum([rec.is_reverse for rec in
                                       self.recList])
        self.ForwardMergedReads = self.MergedReads - self.ReverseMergedReads
        self.ReverseTotalReads = sum([rec.FM for rec in self.recList if
                                      rec.is_reverse])
        self.ForwardTotalReads = self.TotalReads - self.ReverseTotalReads
        try:
            self.AveFamSize = 1. * self.TotalReads / self.MergedReads
        except ZeroDivisionError:
            self.AveFamSize = -1.
        self.TotalAlleleDict = {"A": 0, "C": 0, "G": 0, "T": 0}
        self.SumBQScore = sum(map(oagbq, self.recList))
        self.SumMQScore = sum(map(oagmq, self.recList))
        try:
            self.AveMQ = 1. * self.SumMQScore / lenR
        except ZeroDivisionError:
            self.AveMQ = 0
        try:
            self.AveBQ = 1. * self.SumBQScore / lenR
        except ZeroDivisionError:
            self.AveBQ = 0
        self.consensus = consensus
        self.minMQ = minMQ
        self.minBQ = minBQ
        try:
            self.reverseStrandFraction = (1. * self.ReverseMergedReads /
                                          self.MergedReads)
        except ZeroDivisionError:
            self.reverseStrandFraction = -1.
        try:
            self.MFractionAgreed = nmean(list(map(
                oag("FractionAgreed"),  self.recList)))
        except TypeError:
            pl("Looks like these records have no FractionAgreed attribute. "
               "No worries.", level=logging.DEBUG)
            self.MFractionAgreed = -1.
        self.minFrac = minFracAgreed
        self.minFA = minFA
        try:
            self.MFA = nmean(list(map(oag("FA"), self.recList)))
        except TypeError:
            pl("Looks like these records have no FA attribute. "
               "No worries.", level=logging.DEBUG)
            self.MFA = -1.
        self.FSR = FSR

        # Dealing with transitions (e.g., A-T) and their strandedness
        self.transition = consensus + ">" + self.ALT
        self.strandedTransitionDict = cyfreq(
            ["&&".join([self.transition, is_reverse_to_str(
                rec.is_reverse)]) for rec in self.recList])
        self.StrandCountsDict = {}
        self.StrandCountsDict["reverse"] = sum(map(
            oagir, self.recList))
        self.StrandCountsDict["forward"] = sum(
            rec.is_reverse is False for rec in self.recList)
        self.StrandCountsTotalDict = {}
        self.StrandCountsTotalDict["reverse"] = sum(
            rec.FM for rec in self.recList if rec.is_reverse)
        self.StrandCountsTotalDict["forward"] = sum(
            rec.FM for rec in self.recList if rec.is_reverse is False)
        try:
            self.TotalAlleleFrequency = 1. * self.TotalReads / totalSize
            self.MergedAlleleFrequency = 1. * self.MergedReads / mergedSize
        except ZeroDivisionError:
            self.TotalAlleleFrequency = -1.
            self.MergedAlleleFrequency = -1.
        if(self.ReverseMergedReads != 0 and self.ForwardMergedReads != 0):
            self.BothStrandSupport = True
        else:
            self.BothStrandSupport = False
        if(AABPSD != "default"):
            self.AABPSD = AABPSD
        if(AAMBP != "default"):
            self.AAMBP = AAMBP

        # Check to see if a template supports a variant with both ends
        self.NumberDuplexReads = len([i[0] for i in
                                      cyfreq(map(oag("query_name"),
                                                 self.recList)).iteritems()
                                      if i[1] > 1])
        query_positions = np.array([rec.query_position for
                                    rec in self.recList],
                                   dtype=np.float64)
        self.MBP = nmean(query_positions)
        self.BPSD = nstd(query_positions)
        self.minPVFrac = minPVFrac
        PVFArray = np.array([rec.PVFrac for rec in  self.recList],
                            dtype=np.float64)
        #  PVFArray = [rec.PVFrac for rec in self.recList]

        if(len(PVFArray) == 0):
            self.MPF = -1.
            self.PFSD = -1.
        else:
            self.MPF = nmean(PVFArray)
            self.PFSD = nstd(PVFArray)
        self.maxND = max(rec.opt("ND") for rec in self.recList)


cdef class PCInfo:

    """
    Takes a pysam.PileupColumn covering one base in the reference
    and makes a new class which has "reference" (inferred from
    consensus) and a list of PRData (one for each read).
    The option "duplex required" determines whether or not variants must
    be supported by both reads in a pair for a proper call. Increases
    stringency, might lose some sensitivity for a given sequencing depth.

    exclusionSVTags should be a string of comma-separated SV tags.
    The presence of one of these tags in a read causes it to be thrown out
    of the pileup.

    Note: If a PileupColumn has no reads passing all filters, an error
    AbortMission is thrown so that it can act as a sort of GOTO for rapidly
    getting back to the process without having to remove a None or try passing
    back an invalid object.
    """

    @cython.returns(cystr)
    def getQCFailString(self):
        return ("BQ:%s;FA:%s;MQ:%s" % (self.FailedBQReads,
                                       self.FailedFAReads,
                                       self.FailedMQReads) +
                ";ND:%s;AMP:%s;SV:%s" % (self.FailedNDReads,
                                         self.FailedAMPReads,
                                         self.FailedSVReads) +
                ";AF:%s;QC:%s" % (self.FailedAFReads,
                                  self.FailedQCReads))

    def __init__(self, pPileupColumn_t PileupColumn, int minBQ=0,
                 int minMQ=0, bint requireDuplex=True,
                 float minFracAgreed=0.0, int minFA=0,
                 float minPVFrac=0.66,
                 cystr exclusionSVTags="MDC,LI",
                 bint FracAlignFilter=False, int primerLen=20,
                 cystr experiment="", float minAF=0.25,
                 int maxND=10, object oig1=oig1, object oagir=oagir,
                 object oagqp=oagqp):
        cdef PRInfo_t rec, PRI
        cdef list pileups, fks, svTags, exclusionTagList, discNames
        cdef pPileupRead_t r
        cdef int lenR, rsn
        cdef ndarray[np.float64_t, ndim=1] query_positions
        cdef pPileupRead_t pileupRead
        cdef ndarray[int, ndim=1] SumArray
        cdef dict tmpDict
        pileups = PileupColumn.pileups
        # Get the read pairs which are discordant and get rid of them - one
        # of them has to be wrong!
        if("amplicon" not in experiment):
            primerLen = -1
            # If primerLen < 0, then don't filter by primerLen.
        self.DiscNames = GetDiscReadNames(PileupColumn)
        self.Records = [PRInfo(pileupRead) for
                        pileupRead in pileups if pileupRead.name not in
                        self.DiscNames]
        SumArray = PrunePileupReads(
            self.Records, minMQ=minMQ, minBQ=minBQ, maxND=maxND,
            minFA=minFA, minAF=minAF, primerLen=primerLen)
        self.Records = [PRI for PRI in self.Records if PRI.MQ >= minMQ and PRI.BQ >= minBQ
                        and PRI.FA >= minFA and PRI.AF >= minAF and
                        PRI.query_position >= primerLen and PRI.Pass and PRI.SVPass
                        and PRI.ND > maxND]
        self.FailedMQReads = SumArray[0]
        self.FailedBQReads = SumArray[1]
        self.FailedFAReads = SumArray[2]
        self.FailedAFReads = SumArray[3]
        self.FailedAMPReads = SumArray[4]
        self.FailedSVReads = SumArray[5]
        self.FailedQCReads = SumArray[6]
        self.FailedNDReads = SumArray[7]
        # Remove discordant read pairs.
        self.experiment = experiment
        self.minMQ = minMQ
        self.minBQ = minBQ
        self.contig = PysamToChrInline(PileupColumn.reference_id)
        self.minAF = minAF
        self.pos = PileupColumn.reference_pos
        self.PCol = PileupColumn
        self.excludedSVTagStr = exclusionSVTags
        self.FailedSVReadDict = {}
        self.FailedSVReads = 0

        lenR = len(self.Records)
        if(lenR == 0):
            print(self.getQCFailString())
        rsn = sum(rec.is_reverse for rec in self.Records)
        if(rsn != lenR and rsn != 0):
            self.BothStrandAlignment = True
        else:
            self.BothStrandAlignment = False
        try:
            self.reverseStrandFraction = rsn / (1. * lenR)
        except ZeroDivisionError:
            self.reverseStrandFraction = 0.
        self.MergedReads = lenR
        self.TotalReads = sum([PRI.FM for PRI in self.Records])
        if(lenR > 0):
            self.consensus = sorted(cyfreq(
                map(oagbc, self.Records)).iteritems(),
                key=oig1)[-1][0]
            self.maxND = max(pileupRead.alignment.opt("ND") for
                             pileupRead in pileups)
        else:
            self.consensus = "N"
            pl("Note: Records list empty at position %s" % self.pos)
            self.maxND = 0
            '''
            raise AbortMission("No reads at position passing filters."
                               " Move along - these aren't the "
                               "positions you're looking for.")
            '''
        self.VariantDict = {alt: [rec for rec in self.Records if
                                  rec.BaseCall == alt]
                            for alt in [PRI.BaseCall for PRI in self.Records]}
        query_positions = np.array(map(oagqp, self.Records), dtype=np.float64)
        self.AAMBP = nmean(query_positions)
        self.AABPSD = nstd(query_positions)
        self.AltAlleleData = [AlleleAggregateInfo(
                              self.VariantDict[key],
                              consensus=self.consensus,
                              mergedSize=self.MergedReads,
                              totalSize=self.TotalReads,
                              minMQ=self.minMQ,
                              minBQ=self.minBQ,
                              contig=self.contig,
                              pos=self.pos,
                              DOC=self.MergedReads,
                              DOCTotal=self.TotalReads,
                              NUMALT=len(self.VariantDict),
                              AAMBP=self.AAMBP, AABPSD=self.AABPSD,
                              minFracAgreed=minFracAgreed, minFA=minFA,
                              minPVFrac=minPVFrac, FSR=self.FailedSVReads)
                              for key in self.VariantDict.iterkeys()]
        self.AltAlleleData = [alt for alt in self.AltAlleleData if
                              oag("len")(alt) != 0 and alt is not None]
        # Generates statistics based on transitions, e.g. ref->alt
        self.MergedAlleleDict = {"A": 0, "C": 0, "G": 0, "T": 0}
        self.TotalAlleleDict = {"A": 0, "C": 0, "G": 0, "T": 0}
        for alt in self.AltAlleleData:
            self.MergedAlleleDict[alt.ALT] = alt.MergedReads
            self.TotalAlleleDict[alt.ALT] = alt.TotalReads
        self.TotalCountStr = ",".join(
            [key + ">" + str(self.TotalAlleleDict[key])
             for key in nucList])
        self.MergedCountStr = ",".join(
            [key + ">" + str(self.MergedAlleleDict[key])
             for key in nucList])
        # Process allele frequencies, both by unique and unmerged reads
        self.MergedAlleleFreqDict = {"A": 0., "C": 0., "G": 0., "T": 0.}
        self.TotalAlleleFreqDict = {"A": 0., "C": 0., "G": 0., "T": 0.}
        try:
            for key in self.MergedAlleleDict:
                self.MergedAlleleFreqDict[key] = self.MergedAlleleDict[
                    key] / float(self.MergedReads)
                self.TotalAlleleFreqDict[key] = self.TotalAlleleDict[
                    key] / float(self.TotalReads)
        except ZeroDivisionError:
            for key in self.MergedAlleleDict:
                self.MergedAlleleFreqDict[key] = 0.
                self.TotalAlleleFreqDict[key] = 0.
        self.TotalFracStr = ",".join(
            [key + ">" + str(self.TotalAlleleFreqDict[key])
             for key in nucList])
        self.MergedFracStr = ",".join(
            [key + ">" + str(self.MergedAlleleFreqDict[key])
             for key in nucList])
        # MergedStrandednessRatioDict is the fraction of reverse reads
        # supporting an alternate allele.
        # E.g., if 5 support the alt, but 2 are mapped to the reverse
        # strand, the value is 0.4
        self.MergedStrandednessRatioDict = {"A": 0., "C": 0., "G": 0., "T": 0.}
        self.TotalStrandednessRatioDict = {"A": 0., "C": 0., "G": 0., "T": 0.}
        for alt in self.AltAlleleData:
            self.MergedStrandednessRatioDict[
                alt.ALT] = alt.ReverseMergedReads / alt.MergedReads
            self.TotalStrandednessRatioDict[
                alt.ALT] = alt.ReverseTotalReads / alt.TotalReads
        self.MergedStrandednessStr = "\t".join([
            str(self.MergedStrandednessRatioDict[
                key]) for key in self.MergedStrandednessRatioDict.iterkeys()])
        self.TotalStrandednessStr = "\t".join([
            str(self.TotalStrandednessRatioDict[
                key]) for key in self.TotalStrandednessRatioDict.iterkeys()])
        self.AlleleFreqStr = "\t".join(
            map(str, [self.contig, self.pos, self.consensus,
                      self.MergedReads, self.TotalReads, self.MergedCountStr,
                      self.TotalCountStr, self.MergedFracStr,
                      self.TotalFracStr,
                      self.MergedStrandednessStr,
                      self.TotalStrandednessStr]))

    @cython.returns(AlleleAggregateInfo_t)
    def __getitem__(self, int index):
        return self.AltAlleleData[index]

    @cython.returns(int)
    def __len__(self):
        return len(self.AltAlleleData)

    @cython.returns(str)
    def __str__(self):
        cdef cystr outStr
        outStr = ""
        for alt in self.AltAlleleData:
            outStr += "\t".join(
                np.array([self.contig,
                         self.pos,
                         self.consensus,
                         alt.ALT,
                         alt.TotalReads,
                         alt.MergedReads,
                         alt.TotalAlleleFrequency,
                         alt.MergedAlleleFrequency,
                         alt.StrandCountsTotalDict["reverse"],
                         alt.StrandCountsTotalDict["forward"],
                         alt.StrandCountsDict["reverse"],
                         alt.StrandCountsDict["forward"],
                         self.TotalFracDict[alt.ALT],
                         self.MergedFracDict[alt.ALT],
                         alt.AveFamSize,
                         alt.SumBQScore,
                         alt.AveBQ,
                         alt.SumMQScore,
                         alt.AveMQ]).astype(str)) + "\n"
        self.str = outStr
        return self.str



@cython.returns(tuple)
def PrunepPileupReads(list Records, int minMQ=0, int minBQ=0,
                     int minFA=0, float minAF=0., int maxND=0,
                     int primerLen=-1):
    """
    Returns an array of length 7 - number of failed MQ, BQ, FA, AF,
    amplicon, SV, QC, and ND.
    """
    raise NotImplementedError("Haven't finished updating these objects.")
    cdef py_array retArr = array('i', [0, 0, 0, 0, 0, 0, 0, 0])
    cdef PRInfo_t tpr
    for tpr in Records:
        if tpr.MQ < minMQ:
            retArr[0] += 1
        if tpr.BQ < minBQ:
            retArr[1] += 1
        if tpr.FA < minFA:
            retArr[2] += 1
        if tpr.AF < minAF:
            retArr[3] += 1
        if tpr.query_position < primerLen:
            # Defaults to -1. Gets rid of potentially misprimed nucleotides
            retArr[4] += 1
        if tpr.SVPass is False:
            retArr[5] += 1
        if tpr.Pass is False:
            retArr[6] += 1
        if tpr.ND > maxND:
            retArr[7] += 1
    Records = [tpr for tpr in Records if tpr.MQ >= minMQ and tpr.BQ >= minBQ
               and tpr.FA >= minFA and tpr.AF >= minAF and
               tpr.query_position >= primerLen and tpr.Pass and tpr.SVPass
               and tpr.ND > maxND]
    return retArr, Records


cpdef py_array PrunePileupReads(
        list Records, int minMQ=0, int minBQ=0,
        int minFA=0, float minAF=0., int maxND=20,
        int primerLen=-1):
    """
    Returns an array of length 7 - number of failed MQ, BQ, FA, AF,
    amplicon, SV, QC, and ND.
    """
    cdef ndarray[int, ndim=1] retArr = np.zeros([8], dtype=np.int32)
    cdef PRInfo_t PRI
    for PRI in Records:
        if PRI.MQ < minMQ:
            retArr[0] += 1
        if PRI.BQ < minBQ:
            retArr[1] += 1
        if PRI.FA < minFA:
            retArr[2] += 1
        if PRI.AF < minAF:
            retArr[3] += 1
        if PRI.query_position < primerLen:
            # Defaults to -1. Gets rid of potentially misprimed nucleotides
            retArr[4] += 1
        if PRI.SVPass is False:
            retArr[5] += 1
        if PRI.Pass is False:
            retArr[6] += 1
        if PRI.ND > maxND:
            retArr[7] += 1
    return retArr


class PileupInterval:

    """
    Container for holding summary information over a given interval.
    Written for the pysam PileupColumn data structure.
    Contig should be in the pysam format (IE, a number, not a string)
    """

    def __init__(self, contig="default", start="default",
                 end="default"):
        self.contig = contig
        self.start = int(start)
        self.end = int(end)
        assert self.start < self.end  # Interval must be positive...
        self.TotalCoverage = 0
        self.UniqueCoverage = 0
        self.AvgTotalCoverage = 0
        self.AvgUniqueCoverage = 0

    def updateWithPileupColumn(self, PileupColumn):
        self.end = operator.add(self.end, 1)
        assert self.start < self.end  # Interval must be positive...
        try:
            self.TotalCoverage += sum([int(pu.alignment.opt(
                                      "FM")) for pu in PileupColumn.pileups])
        except KeyError:
            self.TotalCoverage = PileupColumn.nsegments
        self.UniqueCoverage += PileupColumn.nsegments
        self.AvgTotalCoverage = odiv(self.TotalCoverage, float(
            operator.sub(self.end, self.start)))
        self.AvgUniqueCoverage = odiv(self.UniqueCoverage, float(
            operator.sub(self.end, self.start)))

    def __str__(self):
        self.str = "\t".join([str(i) for i in [PysamToChrInline(self.contig),
                                               self.start,
                                               self.end,
                                               self.UniqueCoverage,
                                               self.TotalCoverage,
                                               self.AvgUniqueCoverage,
                                               self.AvgTotalCoverage]])
        return self.str


def BamToCoverageBed(inBAM, outbed="default", mincov=0, minMQ=0, minBQ=0):
    """
    Takes a bam file and creates a bed file containing each position
    """
    pl(("Command required to reproduce this call: "
        "BamToCoverageBed(\"{}\", outbed=".format(inBAM) +
        "\"{}\", mincov={})".format(outbed, mincov)))
    pl(("WARNING: Coverage counts for this script"
        " are wrong. Fix in the works!"
        " It seems to only show up for very long regions."))
    if(outbed == "default"):
        outbed = inBAM[0:-4] + ".doc.bed"
    subprocess.check_call(shlex.split("samtools index {}".format(inBAM)),
                          shell=False)
    if(ospath.isfile(inBAM + ".bai") is False):
        pl("Bam index file was not created. Sorting and indexing.")
        inBAM = HTSUtils.CoorSortAndIndexBam(inBAM)
        pl("Sorted BAM Location: {}".format(inBAM))
    inHandle = pysam.AlignmentFile(inBAM, "rb")
    outHandle = open(outbed, "w")
    workingChr = 0
    workingPos = 0
    outHandle.write("\t".join(["#Chr", "Start", "End",
                               "Number Of Merged Mapped Bases",
                               "Number of Unmerged Mapped Bases",
                               "Avg Merged Coverage",
                               "Avg Total Coverage"]) + "\n")
    pl("Beginning PileupToBed.")
    pileupIterator = inHandle.pileup(max_depth=100000)
    PysamToChrDict = utilBMF.HTSUtils.GetRefIdDicts()["chrtoid"]
    while True:
        try:
            p = pPileupColumn(next(pileupIterator))
        except StopIteration:
            pl("Stopping iterations.")
        PC = PCInfo(p)
        if("Interval" not in locals()):
            if(PC.MergedReads >= mincov):
                Interval = PileupInterval(contig=PC.PCol.reference_id,
                                          start=PC.PCol.reference_pos,
                                          end=operator.add(
                                              PC.PCol.reference_pos, 1))
                try:
                    Interval.updateWithPileupColumn(PC.PCol)
                except AssertionError:
                    del Interval
                workingChr = PC.contig
                workingPos = PC.pos
            continue
        else:
            if(workingChr ==
               PC.PCol.reference_id and PC.PCol.nsegments >= mincov):
                if(operator.sub(PC.PCol.reference_pos, workingPos) == 1):
                    try:
                        Interval.updateWithPileupColumn(PC.PCol)
                        workingPos = operator.add(workingPos, 1)
                    except AssertionError:
                        del Interval
                else:
                    outHandle.write(operator.add(Interval.__str__(), "\n"))
                    if(operator.ge(PC.PCol.nsegments, mincov)):
                        Interval = PileupInterval(
                            contig=PC.PCol.reference_id,
                            start=PC.PCol.reference_pos,
                            end=PC.PCol.reference_pos + 1)
                        workingChr = PC.PCol.reference_id
                        workingPos = PC.PCol.reference_pos
                        Interval.updateWithPileupColumn(PC.PCol)
                    else:
                        del Interval
            else:
                try:
                    Interval.updateWithPileupColumn(PC.PCol)
                    outHandle.write(Interval.__str__() + "\n")
                except AssertionError:
                    del Interval
                if(PC.PCol.nsegments >= mincov):
                    Interval = PileupInterval(contig=PC.PCol.reference_id,
                                              start=PC.PCol.reference_pos,
                                              end=PC.PCol.reference_pos + 1)
                    workingChr = PC.PCol.reference_id
                    workingPos = PC.PCol.reference_pos
                    Interval.updateWithPileupColumn(PC.PCol)
                else:
                    del Interval
    inHandle.close()
    outHandle.close()
    pass

'''
cdef tuple BuildRegionReadDictionaries(AlignmentFile_t handle,
                                       list interval):
    """
    Builds a dictionary from each read in the region for PV and FA.
    """
    cdef AlignedSegment_t tmpRead
    cdef dict PVDict, FADict
    PVDict = {}
    FADict = {}
    for tmpRead in handle.fetch(interval[0], interval[1], interval[2]):
        PVDict[AS_to_key(tmpRead)] =
'''