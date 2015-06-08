# cython: c_string_type=str, c_string_encoding=ascii
# cython: profile=True, cdivision=True, cdivision_warnings=True

from __future__ import division
import subprocess
import os.path
import shlex
import logging
import operator
from operator import attrgetter as oag
from operator import itemgetter as oig
from operator import div as odiv
from itertools import ifilterfalse as iff
try:
    import re2 as re
except ImportError:
    import re

import numpy as np
from numpy import (mean as nmean, max as nmax, sum as nsum,
                   array as nparray, std as nstd)
import pysam
from pysam.calignmentfile import PileupRead as cPileupRead
import cython
from cytoolz import (frequencies as cyfreq,
                     partition as ctpartition)
from itertools import groupby
from utilBMF.ErrorHandling import ThisIsMadness
from utilBMF.HTSUtils import (ReadPair, printlog as pl, pPileupRead,
                              PileupReadPair, ssStringFromRead,
                              PysamToChrDict, nucList, cyStdFlt, cyStdInt)
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


@cython.locals(paired=cython.bint)
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

    def __init__(self, pPileupRead_t PileupRead):
        cdef pysam.calignmentfile.AlignedSegment alignment
        cdef cystr PVSTring
        self.Pass = True
        alignment = PileupRead.alignment
        tags = alignment.tags
        aopt = alignment.opt
        try:
            self.FM = aopt("FM")
        except KeyError:
            self.FM = 1
        """
        except ValueError:
            p = re.compile("\D")
            try:
                self.FM = int(p.sub("", alignment.opt("FM")))
            except ValueError:
                pl("Ain't nothing I can do. Something's wrong"
                   " with your FM tag: {}".format(p.sub(
                       "", alignment.opt("FM"))))
        """
        self.query_name = PileupRead.name
        try:
            self.SVTagString = aopt("SV")
        except KeyError:
            self.SVTagString = "NF"
            # print("Warning: SV Tags unset.")
        self.BaseCall = alignment.query_sequence[
            PileupRead.query_position]
        if(self.BaseCall == "N"):
            self.Pass = False
        self.BQ = alignment.query_qualities[
            PileupRead.query_position]
        self.MQ = alignment.mapq
        self.is_reverse = alignment.is_reverse
        self.is_proper_pair = alignment.is_proper_pair
        self.read = alignment
        self.ssString = ssStringFromRead(alignment)
        self.query_position = PileupRead.query_position
        tagsdictkeys = dict(tags).keys()
        if("FA" in tagsdictkeys):
            self.FA = int(aopt("FA").split(",")[self.query_position])
            self.FractionAgreed = self.FA / (1. * self.FM)
        else:
            self.FA = None
            self.FractionAgreed = None
        self.PV = -1
        self.PV_Array = np.array([])
        self.PVFrac = -1.
        if("PV" in map(oig0, tags)):
            # If there are characters beside digits and commas, then it these
            # values must have been encoded in base 85.
            PVString = aopt("PV")
            self.PV_Array = nparray(PVString.split(','),
                                    dtype=np.int64)
            self.PV = self.PV_Array[self.query_position]
            try:
                self.PVFrac = (1. * self.PV) / nmax(self.PV_Array)
            except ZeroDivisionError:
                pl("ZeroDivision error in PRInfo."
                   "self.PV %s, self.PV_Array %s" % (self.PV, self.PV_Array),
                   level=logging.DEBUG)
                self.PVFrac = -1.
        try:
            self.ND = aopt("ND")
        except KeyError:
            pass
        try:
            self.NF = aopt("NF")
        except KeyError:
            pass
    cpdef object opt(self, cystr arg):
        return self.read.opt(arg)


@cython.returns(cystr)
def is_reverse_to_str(cython.bint boolean):
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
                 object oagmq=oagmq, object nparray=nparray):
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
        if(lenR == 0):
            self = None
            return
        self.len = lenR
        # Total Number of Differences
        if(lenR != 0):
            self.TND = sum(map(oag("ND"), self.recList))
            NFList = nparray(map(oag("NF"), self.recList))
        else:
            self.TND = -1
            NFList = nparray([])
        try:
            self.MNF = nmean(NFList)
            self.maxNF = nmax(NFList)
            self.NFSD = cyStdFlt(NFList)
        except ValueError:
            #  This list must have length zero...
            self.MNF = -1.
            self.maxNF = -1.
            self.NFSD = -1.
        self.TotalReads = sum(map(oag("FM"), self.recList))
        self.MergedReads = lenR
        self.ReverseMergedReads = sum(map(oagir, self.recList))
        self.ForwardMergedReads = self.MergedReads - self.ReverseMergedReads
        self.ReverseTotalReads = sum(map(
            oag("FM"), filter(oagir, self.recList)))
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
        self.ALT = self.recList[0].BaseCall
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
        query_positions = nparray(map(oagqp, self.recList), dtype=np.float64)
        self.MBP = nmean(query_positions)
        self.BPSD = cyStdFlt(query_positions)
        self.minPVFrac = minPVFrac
        PVFArray = nparray(map(oag("PVFrac"), self.recList), dtype=np.float64)
        #  PVFArray = [rec.PVFrac for rec in self.recList]

        if(len(PVFArray) == 0):
            self.MPF = -1.
            self.PFSD = -1.
        else:
            self.MPF = nmean(PVFArray)
            self.PFSD = cyStdFlt(PVFArray)
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
    """

    def __init__(self, pPileupColumn_t PileupColumn, int minBQ=0,
                 int minMQ=0, cython.bint requireDuplex=True,
                 float minFracAgreed=0.0, int minFA=0,
                 float minPVFrac=0.66,
                 cystr exclusionSVTags="MDC,LI",
                 cython.bint FracAlignFilter=False, int primerLen=20,
                 cystr experiment="", float minAF=0.25,
                 int maxND=10, object oig1=oig1, object oagir=oagir,
                 object oagqp=oagqp):
        cdef PRInfo_t rec
        cdef list pileups, fks, svTags, exclusionTagList, discNames
        cdef pPileupRead_t r
        cdef int lenR, rsn
        cdef ndarray[np.float64_t, ndim=1] query_positions
        pileups = PileupColumn.pileups
        # Get the read pairs which are discordant and get rid of them - one
        # of them has to be wrong!
        self.DiscNames = GetDiscReadNames(PileupColumn)
        self.Records = [PRInfo(pileupRead) for
                        pileupRead in pileups
                        if(pileupRead.alignment.mapq >= self.minMQ) and
                        (pileupRead.alignment.query_qualities[
                            pileupRead.query_position] >= self.minBQ) and
                        pileupRead.alignment.opt("FP") == 1 and
                        pileupRead.alignment.opt("FM") >= minFA and
                        pileupRead.alignment.opt("AF") >= minAF and
                        pileupRead.alignment.opt("ND") <= maxND and
                        pileupRead.name not in self.DiscNames]
        # Remove discordant read pairs.
        if("amplicon" in experiment):
            self.ampliconFailed = sum(r for r in pileups
                                      if r.query_position <= primerLen)
            pileups = [r for r in pileups
                       if r.query_position > primerLen]
        self.experiment = experiment
        self.minMQ = minMQ
        #  pl("minMQ: %s" % minMQ)
        self.minBQ = minBQ
        #  pl("minBQ: %s" % minBQ)
        self.contig = PysamToChrDict[PileupColumn.reference_id]
        self.minAF = minAF
        #  pl("Pileup contig: {}".format(self.contig))
        self.pos = PileupColumn.reference_pos
        #  pl("pos: %s" % self.pos)
        self.FailedQCReads = sum(pileupRead.opt("FP") == 0
                                 for pileupRead in pileups)
        self.FailedFMReads = sum(pileupRead.opt("FM") < minFA
                                 for pileupRead in pileups)
        self.FailedAFReads = sum(pileupRead.opt("AF") < minAF
                                 for pileupRead in pileups)
        self.FailedNDReads = sum(pileupRead.opt("ND") > maxND
                                 for pileupRead in pileups)
        self.FailedBQReads = sum(
            pileupRead.alignment.query_qualities[
                pileupRead.query_position] < self.minBQ for
            pileupRead in pileups)
        self.FailedMQReads = sum(
            pileupRead.alignment.mapping_quality < self.minMQ
            for pileupRead in pileups)
        self.PCol = PileupColumn
        self.excludedSVTagStr = exclusionSVTags
        #  pl("Pileup exclusion SV Tags: {}".format(exclusionSVTags))
        svTags = [p.opt("SV") for p in pileups
                  if p.alignment.has_tag("SV")]
        exclusionTagList = exclusionSVTags.split(",")
        svTags = [t for t in svTags if t != "NF" and
                  sum([exTag in t for exTag in exclusionTagList]) != 0]
        self.FailedSVReadDict = {tag: sum([tag in svTag for svTag in svTags])
                                 for tag in exclusionTagList}
        self.FailedSVReads = len(svTags)
        #  pl("Number of reads failed for SV: %s" % self.FailedSVReads)
        """
        if(self.FailedMQReads != 0):
            print("Mapping qualities of reads which failed: {}".format(
                  str([pu.alignment.mapping_quality
                      for pu in pileups
                      if pu.alignment.mapping_quality >= self.minMQ])))
        print("Number of reads failed for BQ < {}: {}".format(self.minBQ,
              self.FailedBQReads))
        print("Number of reads failed for MQ < {}: {}".format(self.minMQ,
              self.FailedMQReads))
        """

        self.Records = filter(oag("Pass"), self.Records)
        lenR = len(self.Records)
        rsn = sum(map(oagir, self.Records))
        if(rsn != lenR and rsn != 0):
            self.BothStrandAlignment = True
        else:
            self.BothStrandAlignment = False
        try:
            self.reverseStrandFraction = sum(
                map(oagir, map(oagread, self.Records))) / (1. * lenR)
        except ZeroDivisionError:
            self.reverseStrandFraction = 0.
        self.MergedReads = lenR
        try:
            self.TotalReads = sum(map(oag("FM"), self.Records))
        except KeyError:
            self.TotalReads = self.MergedReads
        try:
            self.consensus = sorted(cyfreq(
                map(oagbc, self.Records)).iteritems(),
                key=oig1)[-1][0]
        except IndexError:
            self.consensus = "N"  # All bases failed filtered. Oh well.
            pl("Note: PCInfo empty at contig "
               "%s and position %s" % (self.contig, self.pos),
               level=logging.DEBUG)
        self.VariantDict = {alt: [rec for rec in self.Records if
                                  rec.BaseCall == alt]
                            for alt in set(map(oagbc, self.Records))}
        query_positions = nparray(map(oagqp, self.Records), dtype=np.float64)
        self.AAMBP = nmean(query_positions)
        self.AABPSD = cyStdFlt(query_positions)

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
        self.maxND = max(pileupRead.alignment.opt("ND") for
                         pileupRead in pileups)

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
                nparray([self.contig,
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
        self.str = "\t".join([str(i) for i in [PysamToChrDict[self.contig],
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
    if(os.path.isfile(inBAM + ".bai") is False):
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


@cython.locals(MergeDOC=float, TotalDOC=float,
               minMQ=int, minBQ=int)
def CalcWithinBedCoverage(inBAM, bed="default", minMQ=0, minBQ=0,
                          outbed="default"):
    """
    Calculates DOC and creates a bed file containing coverage information
    for each position in a provided bed, only counting reads with
    a given minimum mapping quality or base quality.
    """
    if(outbed == "default"):
        outbed = inBAM[0:-4] + ".doc.bed"
    pl(('Command required to reproduce this call: '
        'CalcWithinBedCoverage("{}", bed='.format(inBAM) +
        '"{}", minMQ="{}", minBQ='.format(bed, minMQ) +
        '"{}", outbed="{}")'.format(minBQ, outbed)))
    if(bed == "default"):
        pl("Bed file required for CalcWithinBedCoverage")
        raise ThisIsMadness("Bedfile required for calculating coverage.")
    bedLines = HTSUtils.ParseBed(bed)
    subprocess.check_call(shlex.split("samtools index {}".format(inBAM)),
                          shell=False)
    if(os.path.isfile(inBAM + ".bai") is False):
        pl("Bam index file could not be created. Sorting and indexing.")
        inBAM = HTSUtils.CoorSortAndIndexBam(inBAM)
        pl("Sorted BAM Location: {}".format(inBAM))
    inHandle = pysam.AlignmentFile(inBAM, "rb")
    outHandle = open(outbed, "w")
    outHandle.write("\t".join(["#Chr", "Start", "End",
                               "Number Of Merged Mapped Bases",
                               "Number of Unmerged Mapped Bases",
                               "Avg Merged Coverage",
                               "Avg Total Coverage",
                               "Mean Family Size",
                               "SD of Family Size"]) + "\n")
    for line in bedLines:
        TotalReads = 0
        MergedReads = 0
        pileupIterator = inHandle.pileup(line[0],
                                         line[1],
                                         line[2],
                                         max_depth=100000,
                                         multiple_iterators=True)
        while True:
            try:
                PC = PCInfo(pPileupColumn(next(pileupIterator)))
            except StopIteration:
                pl("Stopping iteration for bed line: {}".format(line),
                   level=logging.DEBUG)
                break
            TotalReads += PC.TotalReads
            MergedReads += PC.MergedReads
        length = float(operator.sub(line[2], line[1]))
        MergeDOC = odiv(MergedReads, length)
        TotalDOC = odiv(TotalReads, length)
        MeanFamSize = odiv(TotalReads, MergedReads)
        outHandle.write("\t".join(nparray([
            line[0], line[1], line[2],
            MergedReads, TotalReads, MergeDOC,
            TotalDOC]).astype(str)) +
            "\n")
        outHandle.flush()
    inHandle.close()
    outHandle.close()
    return outbed


def CalcWithoutBedCoverage(inBAM, bed="default", minMQ=0, minBQ=0,
                           outbed="default"):
    """
    Calculates DOC and creates a bed file containing each position
    not in provided bed, only counting reads with a given minimum mapping
    quality or base quality.
    """
    if(outbed == "default"):
        outbed = inBAM[0:-4] + ".doc.SBI.bed"
    pl(("Command required to reproduce this call: "
        "CalcWithoutBedCoverage(\"{}\", bed=".format(inBAM) +
        "\"{}\", minMQ=\"{}\", minBQ=".format(bed, minMQ) +
        "\"{}\", outbed={})".format(minBQ, outbed)))
    if(bed == "default"):
        pl("Bed file required for CalcWithoutBedCoverage")
        raise ThisIsMadness("Bedfile required for calculating coverage.")
    bedLines = HTSUtils.ParseBed(bed)
    subprocess.check_call(shlex.split("samtools index {}".format(inBAM)),
                          shell=False)
    if(os.path.isfile(inBAM + ".bai") is False):
        pl("Bam index file could not be created. Sorting and indexing.")
        inBAM = HTSUtils.CoorSortAndIndexBam(inBAM)
        pl("Sorted BAM Location: {}".format(inBAM))
    inHandle = pysam.AlignmentFile(inBAM, "rb")
    outHandle = open(outbed, "w")
    outHandle.write("\t".join(["#Chr", "Start", "End",
                               "Number Of Merged Mapped Bases",
                               "Number of Unmerged Mapped Bases"]) + "\n")
    pileupIterator = inHandle.pileup(max_depth=100000, multiple_iterators=True)
    while True:
        try:
            p = pPileupColumn(next(pileupIterator))
        except StopIteration:
            pl("Finished iterations. (CalcWithoutBedCoverage)")
        PC = PCInfo(p, minMQ=minMQ, minBQ=minBQ)
        if HTSUtils.PosContainedInBed(PC.contig, PC.pos, bedLines):
            continue
        outHandle.write("\t".join(
            [str(i)
             for i in [PC.contig, PC.pos, PC.pos + 1,
                       PC.MergedReads, PC.TotalReads]]) + "\n")
    inHandle.close()
    outHandle.close()
    pass
