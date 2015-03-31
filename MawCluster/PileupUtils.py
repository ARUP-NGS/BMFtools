# cython: c_string_type=str, c_string_encoding=ascii
# cython: profile=True, cdivision=True, cdivision_warnings=True

from __future__ import division
import subprocess
import os.path
import shlex
import logging
import re
import operator
from operator import attrgetter as oag
from operator import itemgetter as oig
from operator import div as odiv
from itertools import ifilterfalse as iff
from collections import Counter

import numpy as np
from numpy import mean as nmean
from numpy import std as nstd
from numpy import max as nmax
from numpy import sum as nsum
from numpy import array as nparray
import pysam
from pysam.calignmentfile import PileupRead as cPileupRead
import cython

from utilBMF.ErrorHandling import ThisIsMadness
from utilBMF.HTSUtils import PysamToChrDict
from utilBMF.HTSUtils import ToStr, ReadPair, printlog as pl
from utilBMF import HTSUtils
import utilBMF


@cython.locals(paired=cython.bint)
def GetDiscordantReadPairs(pPileupColObj):
    """
    Takes a pPileupColumn object (python PileupColumn) as input
    and returns a list of PileupReadPair objects.
    """
    assert isinstance(pPileupColObj, pPileupColumn)
    pileups = pPileupColObj.pileups
    ReadNameCounter = Counter(map(oag("query_name"),
        map(oag("alignment"), pileups)))
    readnames = [i[0] for i in ReadNameCounter.items() if i[1] == 2]
    reads = sorted([read for read in pileups if read.name in readnames],
                   key=lambda x: x.name)
    readpairs = map(PileupReadPair, [reads[2*i:2*i + 2] for i in range(len(reads) // 2)])
    return [pair for pair in readpairs if pair.discordant]


class pPileupRead:
    """
    Python container for the PileupRead proxy in pysam
    """
    def __init__(self, PileupRead):
        self.alignment = PileupRead.alignment
        self.indel = PileupRead.indel
        self.is_del = PileupRead.is_del
        self.level = PileupRead.level
        self.query_position = PileupRead.query_position
        self.name = self.alignment.qname
        self.BaseCall = self.alignment.seq[self.query_position]


class PileupReadPair:

    """
    Holds both bam record objects in a pair of pileup reads.
    Currently, one read unmapped and one read soft-clipped are
    both marked as soft-clipped reads.
    Accepts a list of length two as input.
    """

    def __init__(self, readlist):
        read1, read2 = readlist[0], readlist[1]
        try:
            assert isinstance(read1, cPileupRead) or isinstance(read1, pPileupRead)
        except AssertionError:
            pl("repr(read1): %s" % repr(read1))
            raise ThisIsMadness("PileupReadPair must be initiated with "
                                "pysam.calignmentfile.PileupRead objects or "
                                "pPileup read objects!!")
        try:
            assert len(readlist) == 2
        except AssertionError:
            pl("repr(readlist): %s" % repr(readlist))
            raise ThisIsMadness("readlist must be of length two to make a PileupReadPair!")
        self.RP = ReadPair(read1.alignment, read2.alignment)
        self.read1 = pPileupRead(read1)
        self.read2 = pPileupRead(read2)
        self.discordant = (read1.BaseCall != read2.BaseCall)
        self.name = read1.alignment.query_name
        if(self.discordant):
            if(read1.alignment.is_reverse):
                self.discordanceString = (self.RP.read1_contig + "," +
                                          str(self.read1.alignment.pos -
                                              self.read1.query_position))
            else:
                self.discordanceString = (self.RP.read1_contig + "," +
                                          str(self.read1.alignment.pos +
                                              self.read1.query_position))


class pPileupColumn:
    """
    Python container for the PileupColumn proxy in pysam.
    """
    def __init__(self, PileupColumn):
        self.nsegments = PileupColumn.nsegments
        self.reference_id = PileupColumn.reference_id
        self.reference_pos = PileupColumn.reference_pos
        self.pileups = map(pPileupRead, PileupColumn.pileups)


"""
Contains various utilities for working with pileup generators
and preparing for variant calls.
"""


class PRInfo:

    """
    Created from a pysam.PileupRead object or its python equivalent,
    a pPileupRead.
    Holds family size, SV tags, base quality, mapping quality, and base.
    If any of the "finicky" fields aren't filled (e.g., if BAMs are not
    produced using BMFTools), they are set to None. Check to see if an attribute
    is None first if you'd like to save yourself some headaches.
    """

    def __init__(self, PileupRead):
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
        try:
            self.SVTags = aopt("SV").split(",")
        except KeyError:
            self.SVTags = "NF"
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
        self.ssString = "#".join(
            nparray(sorted(
                [self.read.reference_start,
                 self.read.reference_end])).astype(str))
        self.query_position = PileupRead.query_position
        if("FA" in dict(tags).keys()):
            self.FA_Array = nparray(
                aopt("FA").split(",")).astype(int)
            self.FA = self.FA_Array[self.query_position]
            self.FractionAgreed = odiv(self.FA, float(self.FM))
        else:
            self.FA = None
            self.FractionAgreed = None
            self.FA_Array = None
        p = re.compile("[^0-9,]+")
        self.PV = None
        self.PV_Array = None
        self.PVFrac = None
        if("PV" in map(oig(0), tags)):
            # If there are characters beside digits and commas, then it these
            # values must have been encoded in base 85.
            PVString = aopt("PV")
            try:
                self.PV_Array = nparray(PVString.split(','),
                                         dtype=np.int64)
            except ValueError:
                print("PVString = %s" % PVString)
                raise ValueError("This PV String should only "
                                 "have digits and commas... ???")
            self.PV = self.PV_Array[self.query_position]
            try:
                self.PVFrac = odiv(float(self.PV),
                                           nmax(self.PV_Array))
            except ZeroDivisionError:
                pl("ZeroDivision error in PRInfo."
                   "self.PV %s, self.PV_Array %s" % (self.PV, self.PV_Array))
                self.PVFrac = -1.
        if("ND" in map(oig(0), tags)):
            self.ND = aopt("ND")
        if("NF" in map(oig(0), tags)):
            self.NF = aopt("NF")


def is_reverse_to_str(boolean):
    if(boolean):
        return "reverse"
    elif(boolean is False):
        return "forward"
    else:
        return "unmapped"


class AlleleAggregateInfo:

    """
    Class which holds summary information for a given alt allele
    at a specific position. Meant to be used as part of the PCInfo
    class as values in the VariantDict.
    All alt alleles in this set of reads should be identical.
    recList must be a list of PRInfo objects.

    """

    @cython.locals(minFracAgreed=cython.float, minMQ=cython.long,
                   minBQ=cython.long, minFA=cython.long,
                   minPVFrac=cython.float, FSR=cython.long,
                   lenR=cython.long)
    def __init__(self, recList, consensus="default",
                 mergedSize="default",
                 totalSize="default",
                 minMQ=0,
                 minBQ=0,
                 contig="default",
                 pos="default",
                 DOC="default",
                 DOCTotal="default",
                 NUMALT="default",
                 AABPSD="default", AAMBP="default",
                 minFracAgreed=0.0, minFA=0,
                 minPVFrac=0.5, FSR=-1):
        if(consensus == "default"):
            raise ThisIsMadness("A consensus nucleotide must be provided.")
        if(NUMALT == "default"):
            raise ThisIsMadness("Number of alternate alleles required.")
        else:
            self.NUMALT = NUMALT
        if(DOC == "default"):
            raise ThisIsMadness("Full depth of coverage must be provided.")
        else:
            self.DOC = int(DOC)
        if(DOCTotal != "default"):
            self.DOCTotal = int(DOCTotal)
        if(contig == "default"):
            raise ThisIsMadness("A contig must be provided (string).")
        else:
            self.contig = contig
        if(pos == "default"):
            raise ThisIsMadness("A position (0-based) must be provided (int).")
        else:
            self.pos = int(pos)
        if(mergedSize == "default"):
            raise ThisIsMadness(("mergedSize must be provided: number of "
                                 "PRInfo at given position."))
        if(totalSize == "default"):
            raise ThisIsMadness(("mergedSize must be provided: number of "
                                 "PRInfo at given position."))
        self.recList = [rec for rec in recList
                        if rec.MQ >= minMQ and rec.BQ >= minBQ]
        if("PVFrac" in dir(self.recList[0])):
            self.recList = [rec for rec in self.recList
                            if rec.PVFrac >= minPVFrac
                            or rec.PVFrac is None]
        self.recList = [rec for rec in self.recList
                        if rec.FractionAgreed >= minFracAgreed]

        self.recList = [rec for rec in self.recList if rec.FA >= minFA]
        # Check that all alt alleles are identical
        lenR = len(self.recList)
        self.len = lenR
        # Total Number of Differences
        if(lenR != 0):
            self.TND = sum(map(oag("ND"), self.recList))
            NFList = map(oag("NF"), self.recList)
        else:
            self.TND = -1
            NFList = []
        try:
            self.MNF = nmean(NFList)
            self.maxNF = nmax(NFList)
            self.NFSD = nstd(NFList)
        except ValueError:
            #  This list must have length zero...
            self.MNF = -1.
            self.maxNF = -1.
            self.NFSD = -1.
        try:
            assert(sum([rec.BaseCall == recList[
                0].BaseCall for rec in recList]) == len(recList))
        except AssertionError:
            # print("recList repr: {}".format(repr(recList)))
            raise ThisIsMadness(
                "AlleleAggregateInfo requires that all alt alleles agree.")
        self.TotalReads = nsum(map(oag("FM"), recList))
        self.MergedReads = lenR
        self.ReverseMergedReads = nsum(map(
            oag("is_reverse"), recList))
        self.ForwardMergedReads = self.MergedReads - self.ReverseMergedReads
        self.ReverseTotalReads = nsum(map(
            oag("FM"), filter(oag("is_reverse"), recList)))
        self.ForwardTotalReads = operator.sub(self.TotalReads,
                                              self.ReverseTotalReads)
        try:
            self.AveFamSize = odiv(float(self.TotalReads),
                                           self.MergedReads)
        except ZeroDivisionError:
            self.AveFamSize = -1.
        self.TotalAlleleDict = {"A": 0, "C": 0, "G": 0, "T": 0}
        self.SumBQScore = sum(map(oag("BQ"), recList))
        self.SumMQScore = sum(map(oag("MQ"), recList))
        try:
            self.AveMQ = odiv(float(self.SumMQScore), lenR)
        except ZeroDivisionError:
            self.AveMQ = 0
        try:
            self.AveBQ = odiv(float(self.SumBQScore), lenR)
        except ZeroDivisionError:
            self.AveBQ = 0
        self.ALT = recList[0].BaseCall
        self.consensus = consensus
        self.minMQ = minMQ
        self.minBQ = minBQ
        try:
            self.reverseStrandFraction = odiv(
                float(self.ReverseMergedReads), self.MergedReads)
        except ZeroDivisionError:
            self.reverseStrandFraction = -1.
        try:
            self.MFractionAgreed = nmean(map(
                oag("FractionAgreed"),  recList))
        except TypeError:
            pl("Looks like these records have no FractionAgreed attribute. "
               "No worries.", level=logging.DEBUG)
            self.MFractionAgreed = -1.
        self.minFrac = minFracAgreed
        self.minFA = minFA
        try:
            self.MFA = nmean(map(oag("FA"), recList))
        except TypeError:
            pl("Looks like these records have no FA attribute. "
               "No worries.", level=logging.DEBUG)
            self.MFA = -1.
        self.FSR = FSR

        # Dealing with transitions (e.g., A-T) and their strandedness
        self.transition = ">".join([consensus, self.ALT])
        self.strandedTransitions = {}
        self.strandedTransitionDict = Counter(
            ["&&".join([self.transition, is_reverse_to_str(
                rec.is_reverse)]) for rec in self.recList])
        self.strandedTotalTransitionDict = Counter(
            ["&&".join([self.transition, is_reverse_to_str(
                rec.is_reverse)]) for rec in self.recList] * rec.FM)
        self.strandedMergedTransitionDict = Counter(
            ["&&".join([self.transition, is_reverse_to_str(
                rec.is_reverse)]) for rec in self.recList])
        self.StrandCountsDict = {}
        self.StrandCountsDict["reverse"] = sum(map(
            oag("is_reverse"), self.recList))
        self.StrandCountsDict["forward"] = sum([
            rec.is_reverse is False for rec in self.recList])
        self.StrandCountsTotalDict = {}
        self.StrandCountsTotalDict["reverse"] = sum(
            [rec.FM for rec in self.recList if rec.is_reverse])
        self.StrandCountsTotalDict["forward"] = sum(
            [rec.FM for rec in self.recList if rec.is_reverse is False])
        try:
            self.TotalAlleleFrequency = self.TotalReads / totalSize
            self.MergedAlleleFrequency = self.MergedReads / mergedSize
        except ZeroDivisionError:
            self.TotalAlleleFrequency = -1.0
            self.MergedAlleleFrequency = -1.0
        if(self.ReverseMergedReads != 0 and self.ForwardMergedReads != 0):
            self.BothStrandSupport = True
        else:
            self.BothStrandSupport = False
        if(AABPSD != "default"):
            self.AABPSD = AABPSD
        if(AAMBP != "default"):
            self.AAMBP = AAMBP

        # Check to see if a read pair supports a variant with both ends
        ReadNameCounter = Counter(map(
            oag("query_name"),
            map(oag("read"), recList)))
        self.NumberDuplexReads = sum([
            ReadNameCounter[key] > 1 for key in ReadNameCounter.keys()])
        query_positions = nparray(map(
            oag("query_position"), recList)).astype(float)
        self.MBP = nmean(query_positions)
        self.BPSD = nstd(query_positions)
        self.minPVFrac = minPVFrac
        PVFArray = [i.PVFrac for i in self.recList if i.PVFrac is not None]
        if(len(PVFArray) == 0):
            self.MPF = -1.
            self.PFSD = -1.
        else:
            self.MPF = nmean(PVFArray)
            self.PFSD = nstd(PVFArray)


class PCInfo:

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

    @cython.locals(reverseStrandFraction=cython.float,
                   BothStrandAlignment=cython.bint,
                   primerLen=cython.long, minBQ=cython.long,
                   minMQ=cython.long, minPVFrac=cython.float)
    def __init__(self, PileupColumn, minBQ=0, minMQ=0,
                 requireDuplex=True,
                 minFracAgreed=0.0, minFA=0, minPVFrac=0.66,
                 exclusionSVTags="MDC,LI,ORB",
                 FracAlignFilter=False, primerLen=20,
                 experiment=""):
        assert isinstance(PileupColumn, pPileupColumn)
        pileups = PileupColumn.pileups
        if("amplicon" in experiment):
            self.ampliconFailed = sum([r for r in pileups
                                       if r.query_position <= primerLen])
            pileups = [r for r in pileups
                       if r.query_position > primerLen]
        self.experiment = experiment
        self.minMQ = minMQ
        #  pl("minMQ: %s" % minMQ)
        self.minBQ = minBQ
        #  pl("minBQ: %s" % minBQ)
        from collections import Counter
        self.contig = PysamToChrDict[PileupColumn.reference_id]
        #  pl("Pileup contig: {}".format(self.contig))
        self.pos = PileupColumn.reference_pos
        #  pl("pos: %s" % self.pos)
        self.FailedQCReads = sum([pileupRead.alignment.opt("FP") == 0
                                  for pileupRead in pileups])
        self.FailedFMReads = sum([pileupRead.alignment.opt("FM") < minFA
                                  for pileupRead in pileups])
        self.FailedBQReads = sum(
            [pileupRead.alignment.query_qualities
             [pileupRead.query_position] < self.minBQ
             for pileupRead in pileups])
        self.FailedMQReads = sum(
            [pileupRead.alignment.mapping_quality < self.minMQ
             for pileupRead in pileups])
        self.PCol = PileupColumn
        self.excludedSVTagStr = exclusionSVTags
        #  pl("Pileup exclusion SV Tags: {}".format(exclusionSVTags))
        self.FailedSVReadDict = {tag: sum([p.alignment.has_tag("SV") and tag
                                           not in p.alignment.opt("SV") for
                                           p in pileups])
                                 for tag in exclusionSVTags.split(",")}
        """
        self.FailedSVReadDict = {}
        for tag in self.excludedSVTags:
            self.FailedSVReadDict[tag] = sum([p.alignment.has_tag("SV") and tag
                                              not in p.alignment.opt("SV")
                                              for p in pileups])
        """
        self.FailedSVReads = sum([self.FailedSVReadDict[key] for key
                                  in self.FailedSVReadDict.keys()])
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
        self.Records = [PRInfo(pileupRead) for pileupRead in pileups
                        if(pileupRead.alignment.mapq >= self.minMQ) and
                        (pileupRead.alignment.query_qualities[
                            pileupRead.query_position] >= self.minBQ) and
                        pileupRead.alignment.opt("FP") == 1 and
                        pileupRead.alignment.opt("FM") >= minFA]
        self.Records = filter(oag("Pass"), self.Records)
        lenR = len(self.Records)
        rsn = sum(map(oag("is_reverse"), self.Records))
        if(rsn != lenR and rsn != 0):
            self.BothStrandAlignment = True
        else:
            self.BothStrandAlignment = False
        try:
            self.reverseStrandFraction = sum(map(
                oag("is_reverse"), map(oag("read"),
                                       self.Records))) / float(lenR)
        except ZeroDivisionError:
            self.reverseStrandFraction = 0.
        self.MergedReads = lenR
        try:
            self.TotalReads = sum(map(oag("FM"), self.Records))
        except KeyError:
            self.TotalReads = self.MergedReads
        try:
            self.consensus = Counter(
                map(oag("BaseCall"), self.Records)).most_common(1)[0][0]
        except IndexError:
            self.consensus = Counter(map(
                oag("BaseCall"),
                map(PRInfo, pileups))).most_common(1)[0][0]
        self.VariantDict = {}
        for alt in list(set(map(oag("BaseCall"), self.Records))):
            self.VariantDict[alt] = [
                rec for rec in self.Records if rec.BaseCall == alt]
        query_positions = map(oag("query_position"),
                              self.Records)
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
                              for key in self.VariantDict.keys()]
        self.AltAlleleData = [i for i in self.AltAlleleData if
                              oag("len")(i) != 0]
        self.TotalFracDict = {"A": 0., "C": 0., "G": 0., "T": 0.}
        for alt in self.AltAlleleData:
            self.TotalFracDict[
                alt.ALT] = 1. * alt.TotalReads / self.TotalReads
        self.TotalFracStr = ",".join(
            [">".join([key, str(self.TotalFracDict[key])])
             for key in self.TotalFracDict.keys()])
        self.TotalCountDict = {"A": 0, "C": 0, "G": 0, "T": 0}
        for alt in self.AltAlleleData:
            self.TotalCountDict[
                alt.ALT] = alt.TotalReads
        self.TotalCountStr = ",".join(
            [">".join([key, str(self.TotalCountDict[key])])
             for key in self.TotalCountDict.keys()])
        self.MergedFracDict = {"A": 0., "C": 0., "G": 0., "T": 0.}
        for alt in self.AltAlleleData:
            self.MergedFracDict[
                alt.ALT] = odiv(float(alt.MergedReads),
                                        self.MergedReads)
        self.MergedFracStr = ",".join(
            [">".join([key, str(self.MergedFracDict[key])])
             for key in self.MergedFracDict.keys()])
        self.MergedCountDict = {"A": 0, "C": 0, "G": 0, "T": 0}
        for alt in self.AltAlleleData:
            self.MergedCountDict[
                alt.ALT] = alt.MergedReads
        self.MergedCountStr = ",".join(
            [">".join([key, str(self.MergedCountDict[key])])
             for key in self.MergedCountDict.keys()])
        # Generates statistics based on transitions, e.g. ref->alt
        TransMergedCounts = {}
        for alt in self.AltAlleleData:
            try:
                TransMergedCounts[alt.transition] = operator.add(
                    TransMergedCounts[alt.transition], alt.MergedReads)
            except KeyError:
                TransMergedCounts[alt.transition] = alt.MergedReads
        self.TransMergedCounts = TransMergedCounts
        TransTotalCounts = {}
        for alt in self.AltAlleleData:
            try:
                TransTotalCounts[alt.transition] = operator.add(
                    alt.TotalReads, TransTotalCounts[alt.transition])
            except KeyError:
                TransTotalCounts[alt.transition] = alt.TotalReads
        self.TransTotalCounts = TransTotalCounts
        self.StrandedTransTotalCounts = {}
        for alt in self.AltAlleleData:
            for trans in alt.strandedTotalTransitionDict.keys():
                try:
                    self.StrandedTransTotalCounts[
                        trans] += alt.strandedTotalTransitionDict[trans]
                except KeyError:
                    self.StrandedTransTotalCounts[
                        trans] = alt.strandedTotalTransitionDict[trans]
        self.StrandedTransMergedCounts = {}
        for alt in self.AltAlleleData:
            for trans in alt.strandedMergedTransitionDict.keys():
                try:
                    self.StrandedTransMergedCounts[
                        trans] += alt.strandedMergedTransitionDict[trans]
                except KeyError:
                    self.StrandedTransMergedCounts[
                        trans] = alt.strandedMergedTransitionDict[trans]
        self.MergedAlleleDict = {"A": 0, "C": 0, "G": 0, "T": 0}
        self.TotalAlleleDict = {"A": 0, "C": 0, "G": 0, "T": 0}
        for alt in self.AltAlleleData:
            self.MergedAlleleDict[alt.ALT] = alt.MergedReads
            self.TotalAlleleDict[alt.ALT] = alt.TotalReads
        # Process allele frequencies, both by unique and unmerged reads
        self.MergedAlleleFreqDict = {"A": 0., "C": 0., "G": 0., "T": 0.}
        self.TotalAlleleFreqDict = {"A": 0., "C": 0., "G": 0., "T": 0.}
        try:
            for key in self.MergedAlleleDict:
                self.MergedAlleleFreqDict[key] = odiv(
                    self.MergedAlleleDict[key],
                    float(self.MergedReads))
                self.TotalAlleleFreqDict[key] = odiv(
                    self.TotalAlleleDict[key],
                    float(self.TotalReads))
        except ZeroDivisionError:
            for key in self.MergedAlleleDict:
                self.MergedAlleleFreqDict[key] = 0.
                self.TotalAlleleFreqDict[key] = 0.
        self.MergedAlleleCountStr = "\t".join(
            map(ToStr, [self.MergedAlleleDict["A"],
                        self.MergedAlleleDict["C"],
                        self.MergedAlleleDict["G"],
                        self.MergedAlleleDict["T"]]))
        self.TotalAlleleCountStr = "\t".join(
            map(ToStr, [self.TotalAlleleDict["A"],
                        self.TotalAlleleDict["C"],
                        self.TotalAlleleDict["G"],
                        self.TotalAlleleDict["T"]]))
        self.MergedAlleleFreqStr = "\t".join(
            map(ToStr, [self.MergedAlleleFreqDict["A"],
                        self.MergedAlleleFreqDict["C"],
                        self.MergedAlleleFreqDict["G"],
                        self.MergedAlleleFreqDict["T"]]))
        self.TotalAlleleFreqStr = "\t".join(
            map(ToStr, [self.TotalAlleleFreqDict["A"],
                        self.TotalAlleleFreqDict["C"],
                        self.TotalAlleleFreqDict["G"],
                        self.TotalAlleleFreqDict["T"]]))
        # MergedStrandednessRatioDict is the fraction of reverse reads
        # supporting an alternate allele.
        # E.g., if 5 support the alt, but 2 are mapped to the reverse
        # strand, the value is 0.4
        self.MergedStrandednessRatioDict = {"A": 0., "C": 0., "G": 0., "T": 0.}
        self.TotalStrandednessRatioDict = {"A": 0., "C": 0., "G": 0., "T": 0.}
        for alt in self.AltAlleleData:
            self.MergedStrandednessRatioDict[
                alt.ALT] = odiv(alt.ReverseMergedReads,
                                        float(alt.MergedReads))
            self.TotalStrandednessRatioDict[
                alt.ALT] = odiv(alt.ReverseTotalReads,
                                        float(alt.TotalReads))
        self.MergedStrandednessStr = "\t".join([
            str(self.MergedStrandednessRatioDict[
                key]) for key in self.MergedStrandednessRatioDict.keys()])
        self.TotalStrandednessStr = "\t".join([
            str(self.TotalStrandednessRatioDict[
                key]) for key in self.TotalStrandednessRatioDict.keys()])
        self.AlleleFreqStr = "\t".join(
            [str(i) for i in [self.contig,
                              self.pos,
                              self.consensus,
                              self.MergedReads,
                              self.TotalReads,
                              self.MergedAlleleCountStr,
                              self.TotalAlleleCountStr,
                              self.MergedAlleleFreqStr,
                              self.TotalAlleleFreqStr,
                              self.MergedStrandednessStr,
                              self.TotalStrandednessStr]])

    def ToString(self, header=False):
        outStr = ""
        if(header):
            outStr = ("#Chr\tPos (0-based)\tRef (Consensus)\tAlt\tTotal "
                      "Reads\tMerged Reads\tTotal Allele Frequency\tMerged "
                      "Allele Frequency\tReverse Total Reads\tForward Total"
                      " Reads\tReverse Merged Reads\tForward Merged Reads"
                      "\tFraction Of Total Reads\t"
                      "Fraction Of Merged Reads\tAverage "
                      "Family Size\t"
                      "BQ Sum\tBQ Mean\tMQ Sum\tMQ Mean\n")
        for alt in self.AltAlleleData:
            outStr += "\t".join(
                map(ToStr, [self.contig,
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
                            alt.AveMQ])) + "\n"
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
        self.str = self.ToString()

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

    def ToString(self):
        self.str = "\t".join([str(i) for i in [PysamToChrDict[self.contig],
                                               self.start,
                                               self.end,
                                               self.UniqueCoverage,
                                               self.TotalCoverage,
                                               self.AvgUniqueCoverage,
                                               self.AvgTotalCoverage]])
        return self.str


def CustomPileupFullGenome(inputBAM,
                           PileupTsv="default",
                           TransitionTable="default",
                           StrandedTTable="default",
                           progRepInterval=10000,
                           minBQ=0,
                           minMQ=0):
    """
    A pileup tool for creating a tsv for each position in the bed file.
    Used for calling SNPs with high confidence.
    Also creates several tables:
    1. Counts for all consensus->alt transitions (By Total and Merged reads)
    2. Counts for the above, specifying strandedness
    3. Number of Merged Reads supporting each allele
    """
    TransTotalDict = {}
    TransMergedDict = {}
    StrandedTransTotalDict = {}
    StrandedTransMergedDict = {}
    NumTransitionsTotal = 0
    NTransMerged = 0
    MergedReadsProcessed = 0
    TotalReadsProcessed = 0
    if(PileupTsv == "default"):
        PileupTsv = inputBAM[0:-4] + ".Pileup.tsv"
    if(TransitionTable == "default"):
        TransitionTable = inputBAM[0:-4] + ".Trans.tsv"
    if(StrandedTTable == "default"):
        StrandedTTable = inputBAM[0:-4] + ".StrandedTrans.tsv"
    subprocess.check_call(["samtools", "index", inputBAM])
    if(os.path.isfile(inputBAM + ".bai") is False):
        pl("Couldn\'t index BAM - coor sorting, then indexing!")
        inputBAM = HTSUtils.CoorSortAndIndexBam(inputBAM, uuid=True)
    NumProcessed = 0  # Number of processed positions in pileup
    bamHandle = pysam.AlignmentFile(inputBAM, "rb")
    PileupHandle = open(PileupTsv, "w")
    TransHandle = open(TransitionTable, "w")
    StrandedTransHandle = open(StrandedTTable, "w")
    FirstLine = True
    pileupIterator = bamHandle.pileup(max_depth=100000, multiple_iterators=True)
    p = pPileupColumn(next(pileupIterator))
    for pileupColumn in p.pileups:
        NumProcessed = operator.add(NumProcessed, 1)
        if((NumProcessed) % progRepInterval == 0):
            pl("Number of positions processed: {}".format(
                NumProcessed))
            pl("Total reads processed: {}".format(TotalReadsProcessed))
            pl("Merged reads processed: {}".format(MergedReadsProcessed))
        PColSum = PCInfo(pileupColumn, minBQ=minBQ, minMQ=minMQ)
        MergedReadsProcessed = operator.add(MergedReadsProcessed,
                                            PColSum.MergedReads)
        TotalReadsProcessed = operator.add(TotalReadsProcessed,
                                           PColSum.TotalReads)
        if(FirstLine):
            PileupHandle.write(PColSum.ToString(header=True))
            FirstLine = False
        else:
            PileupHandle.write(PColSum.ToString())
        for key in PColSum.TransMergedCounts.keys():
            try:
                TransMergedDict[
                    key] += PColSum.TransMergedCounts[key]
            except KeyError:
                TransMergedDict[
                    key] = PColSum.TransMergedCounts[key]
            NTransMerged += PColSum.TransMergedCounts[
                key]
        for key in PColSum.TransTotalCounts.keys():
            try:
                TransTotalDict[
                    key] += PColSum.TransTotalCounts[key]
            except KeyError:
                TransTotalDict[
                    key] = PColSum.TransTotalCounts[key]
            NumTransitionsTotal += PColSum.TransTotalCounts[
                key]
        for key in PColSum.StrandedTransMergedCounts.keys():
            try:
                StrandedTransMergedDict[
                    key] += PColSum.StrandedTransMergedCounts[
                        key]
            except KeyError:
                StrandedTransMergedDict[
                    key] = PColSum.StrandedTransMergedCounts[
                        key]
        for key in PColSum.StrandedTransTotalCounts.keys():
            try:
                StrandedTransTotalDict[
                    key] += PColSum.StrandedTransTotalCounts[
                        key]
            except KeyError:
                StrandedTransTotalDict[
                    key] = PColSum.StrandedTransTotalCounts[
                        key]
    TransHandle.write(("Transition\tTotal Reads With Transition (Unflattened)"
                       "\tMerged Reads With Transition\tFraction Of Total "
                       "Transitions\tFraction Of Merged Transitions\n"))
    for key in TransTotalDict.keys():
        TransHandle.write("{}\t{}\t{}\n".format(key,
                                                TransTotalDict[key],
                                                TransMergedDict[key],
                                                odiv(
                                                    TransTotalDict[key],
                                                    float(
                                                        NumTransitionsTotal)),
                                                odiv(
                                                    TransMergedDict[key],
                                                    float(
                                                        NTransMerged))))
    StrandedTransHandle.write(("Transition+Strandedness\tTotal Reads "
                               "(Unflattened)\tMergedReads With Transition\t"
                               "Fraction Of Total (Unflattened) Transitions"
                               "\tFraction of Merged Transitions\n"))
    for key in StrandedTransTotalDict.keys():
        StrandedTransHandle.write(
            "{}\t{}\t{}\n".format(
                key,
                StrandedTransTotalDict[key],
                StrandedTransMergedDict[key],
                odiv(StrandedTransTotalDict[key],
                             float(NumTransitionsTotal)),
                odiv(StrandedTransMergedDict[key],
                             float(NTransMerged)),
            ))
    pl("Transition Table: {}".format(TransitionTable))
    pl("Stranded Transition Table: {}".format(StrandedTTable))
    TransHandle.close()
    PileupHandle.close()
    return PileupTsv


def CustomPileupToTsv(inputBAM,
                      PileupTsv="default",
                      TransitionTable="default",
                      StrandedTTable="default",
                      bedfile="default",
                      progRepInterval=1000,
                      CalcAlleleFreq=True,
                      minMQ=0,
                      minBQ=0):
    """
    A pileup tool for creating a tsv for each position in the bed file.
    Used for calling SNPs with high confidence.
    Also creates several tables:
    1. Counts for all consensus->alt transitions (By Total and Merged reads)
    2. Counts for the above, specifying strandedness
    3. Number of Merged Reads supporting each allele
    """
    TransTotalDict = {}
    TransMergedDict = {}
    StrandedTransTotalDict = {}
    StrandedTransMergedDict = {}
    NumTransitionsTotal = 0
    NTransMerged = 0
    MergedReadsProcessed = 0
    TotalReadsProcessed = 0
    if(PileupTsv == "default"):
        PileupTsv = inputBAM[0:-4] + ".Pileup.tsv"
    if(TransitionTable == "default"):
        TransitionTable = inputBAM[0:-4] + ".Trans.tsv"
    if(StrandedTTable == "default"):
        StrandedTTable = inputBAM[0:-4] + ".StrandedTrans.tsv"
    if(bedfile == "default"):
        return CustomPileupFullGenome(inputBAM, PileupTsv=PileupTsv,
                                      TransitionTable=TransitionTable,
                                      StrandedTTable=StrandedTTable,
                                      progRepInterval=progRepInterval)
    bedlines = HTSUtils.ParseBed(bedfile)
    try:
        subprocess.check_call(["samtools", "index", inputBAM])
    except subprocess.CalledProcessError:
        pl("Couldn't index BAM - coor sorting, then indexing!")
        inputBAM = HTSUtils.CoorSortAndIndexBam(inputBAM, uuid=True)
    NumPos = sum([operator.sub(line[2], line[1]) for line in bedlines])
    NumProcessed = 0  # Number of processed positions in pileup
    pl("Number of positions in bed file: {}".format(NumPos))
    bamHandle = pysam.AlignmentFile(inputBAM, "rb")
    PileupHandle = open(PileupTsv, "w")
    TransHandle = open(TransitionTable, "w")
    StrandedTransHandle = open(StrandedTTable, "w")
    FirstLine = True
    for line in bedlines:
        pileupIterator = bamHandle.pileup(line[0], line[1], line[2],
                                          max_depth=100000)
        p = pPileupColumn(next(pileupIterator))
        for pileupColumn in p.pileups:
            NumProcessed += 1
            if((NumProcessed) % progRepInterval == 0):
                pl("Number of positions processed: {}".format(
                    NumProcessed + 1))
                pl("{:.1%} complete".format(NumProcessed / float(NumPos)))
                pl("Total reads processed: {}".format(TotalReadsProcessed))
                pl("Merged reads processed: {}".format(
                    MergedReadsProcessed))
            PColSum = PCInfo(pileupColumn, minMQ=minMQ, minBQ=minBQ)
            MergedReadsProcessed += PColSum.MergedReads
            TotalReadsProcessed += PColSum.TotalReads
            if(FirstLine):
                PileupHandle.write(PColSum.ToString(header=True))
                FirstLine = False
            else:
                PileupHandle.write(PColSum.ToString())
            for key in PColSum.TransMergedCounts.keys():
                try:
                    TransMergedDict[
                        key] += PColSum.TransMergedCounts[key]
                except KeyError:
                    TransMergedDict[
                        key] = PColSum.TransMergedCounts[key]
                NTransMerged += PColSum.TransMergedCounts[
                    key]
            for key in PColSum.TransTotalCounts.keys():
                try:
                    TransTotalDict[
                        key] += PColSum.TransTotalCounts[key]
                except KeyError:
                    TransTotalDict[
                        key] = PColSum.TransTotalCounts[key]
                NumTransitionsTotal += PColSum.TransTotalCounts[
                    key]
            for key in PColSum.StrandedTransMergedCounts.keys():
                try:
                    StrandedTransMergedDict[
                        key] += PColSum.StrandedTransMergedCounts[
                            key]
                except KeyError:
                    StrandedTransMergedDict[
                        key] = PColSum.StrandedTransMergedCounts[
                            key]
            for key in PColSum.StrandedTransTotalCounts.keys():
                try:
                    StrandedTransTotalDict[
                        key] += PColSum.StrandedTransTotalCounts[
                            key]
                except KeyError:
                    StrandedTransTotalDict[
                        key] = PColSum.StrandedTransTotalCounts[
                            key]
    TransHandle.write((
        "Transition\tTotal Reads With Transition (Unflatt"
        "ened)\tMerged Reads With Transition\tFraction Of Total "
        "Transitions\tFraction Of Merged Transitions\n"))
    for key in TransTotalDict.keys():
        if(key[0] != key[3]):
            TransHandle.write(
                "{}\t{}\t{}\n".format(
                    key,
                    TransTotalDict[key],
                    TransMergedDict[key],
                    odiv(TransTotalDict[key],
                                 float(NumTransitionsTotal)),
                    odiv(TransMergedDict[key],
                                 float(NTransMerged))))
    StrandedTransHandle.write(("Transition+Strandedness\tTotal Reads "
                               "(Unflattened)\tMergedReads With Transition\t"
                               "Fraction Of Total (Unflattened) Transitions"
                               "\tFraction of Merged Transitions\n"))
    for key in StrandedTransTotalDict.keys():
        if(key.split("&&")[0].split(">")[0] != key.split(
                "&&")[0].split(">")[1]):
            StrandedTransHandle.write(
                "{}\t{}\t{}\t{}\t{}\n".format(
                    key,
                    StrandedTransTotalDict[key],
                    StrandedTransMergedDict[key],
                    StrandedTransTotalDict[key] /
                    float(NumTransitionsTotal),
                    StrandedTransMergedDict[key] /
                    float(NTransMerged),
                ))
    pl("Transition Table: {}".format(TransitionTable))
    pl("Stranded Transition Table: {}".format(StrandedTTable))
    TransHandle.close()
    PileupHandle.close()
    if(CalcAlleleFreq):
        AlleleFreqTsv = AlleleFrequenciesByBase(inputBAM,
                                                bedfile=bedfile,
                                                minMQ=minMQ,
                                                minBQ=minBQ)
        pl("Optional allele frequency table generated. Path: {}".format(
            AlleleFreqTsv))
    return PileupTsv


def AlleleFrequenciesByBase(inputBAM,
                            outputTsv="default",
                            progRepInterval=10000,
                            minMQ=0,
                            minBQ=0,
                            bedfile="default"):
    """
    Creates a tsv file with counts and frequencies for each allele at
    each position. I should expand this to include strandedness information.
    """
    pl(("Command required to reproduce results: "
        "AlleleFrequenciesByBase(inputBAM=\"{}\",".format(inputBAM) +
        " outputTsv=\"{}\", progRepInterval=".format(outputTsv) +
        "\"{}\", minMQ=\"{}\", )".format(progRepInterval, minMQ) +
        "minBQ=\"{}\"".format(minBQ)))
    if(outputTsv == "default"):
        outputTsv = inputBAM[0:-4] + ".allele.freq.tsv"
    try:
        subprocess.check_call(["samtools", "index", inputBAM])
    except subprocess.CalledProcessError:
        pl("Couldn't index BAM - coor sorting, then indexing!")
        inputBAM = HTSUtils.CoorSortAndIndexBam(inputBAM, uuid=True)
    NumProcessed = 0  # Number of processed positions in pileup
    inHandle = pysam.AlignmentFile(inputBAM, "rb")
    outHandle = open(outputTsv, "w")
    outHandle.write("\t".join(["#Contig",
                               "Position (0-based)",
                               "Consensus Base",
                               "Merged DOC",
                               "Total DOC",
                               "Merged Count: A",
                               "Merged Count: C",
                               "Merged Count: G",
                               "Merged Count: T",
                               "Total Count: A",
                               "Total Count: C",
                               "Total Count: G",
                               "Total Count: T",
                               "Merged Freq: A",
                               "Merged Freq: C",
                               "Merged Freq: G",
                               "Merged Freq: T",
                               "Total Freq: A",
                               "Total Freq: C",
                               "Total Freq: G",
                               "Total Freq: T",
                               "Merged Reverse fraction: A",
                               "Merged Reverse fraction: C",
                               "Merged Reverse fraction: G",
                               "Merged Reverse fraction: T",
                               "Total Reverse fraction: A",
                               "Total Reverse fraction: C",
                               "Total Reverse fraction: G",
                               "Total Reverse fraction: T",
                               ]) + "\n")
    if(bedfile == "default"):
        pileupIterator = inHandle.pileup(max_depth=100000)
        p = pPileupColumn(next(pileupIterator))
        for pileup in p.pileups:
            NumProcessed += 1
            if(NumProcessed % progRepInterval == 0):
                pl("Number of base positions processed: {}".format(
                    NumProcessed))
            PColInfo = PCInfo(pileup, minMQ=minMQ, minBQ=minBQ)
            outHandle.write(PColInfo.AlleleFreqStr + "\n")
    else:
        bedLines = HTSUtils.ParseBed(bedfile)
        for line in bedLines:
            # print("Now running through the bedLine: {}".format(line))
            puIterator = inHandle.pileup(reference=line[0], start=line[1],
                                         end=line[2],
                                         max_depth=100000)
            p = pPileupColumn(next(puIterator))
            for pileup in p.pileups:
                NumProcessed += 1
                if(NumProcessed % progRepInterval == 0):
                    pl("Number of base positions processed: {}".format(
                        NumProcessed))
                outHandle.write(PCInfo(pileup,
                                       minMQ=minMQ,
                                       minBQ=minBQ).AlleleFreqStr + "\n")
    inHandle.close()
    outHandle.close()
    return outputTsv


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
    ChrToPysamDict = utilBMF.HTSUtils.GetRefIdDicts()["chrtoid"]
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
                    outHandle.write(operator.add(Interval.ToString(), "\n"))
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
                    outHandle.write(Interval.ToString() + "\n")
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


@cython.locals(MergeDOC=cython.float, TotalDOC=cython.float,
               minMQ=cython.long, minBQ=cython.long)
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
        outHandle.write(
            "\t".join(
                  nparray([line[0],
                  line[1],
                  line[2],
                  MergedReads, TotalReads, MergeDOC,
                  TotalDOC,
                            ]).astype(str))
                        + "\n")
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
