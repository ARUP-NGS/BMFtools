import subprocess
import numpy as np
import os.path
import shlex

import pysam

from utilBMF.HTSUtils import ThisIsMadness, printlog as pl, PysamToChrDict
from utilBMF import HTSUtils
import utilBMF

"""
Contains various utilities for working with pileup generators
and preparing for variant calls.
"""


class PRInfo:

    """
    Created from a pysam.PileupRead object.
    Holds family size, SV tags, base quality,
    mapping quality, and base.
    """

    def __init__(self, PileupRead):
        try:
            self.FM = int(PileupRead.alignment.opt("FM"))
        except KeyError:
            self.FM = 1
        try:
            self.SVTags = PileupRead.alignment.opt("SV").split(",")
        except KeyError:
            self.SVTags = "NF"
            # print("Warning: SV Tags unset.")
        self.BaseCall = PileupRead.alignment.query_sequence[
            PileupRead.query_position]
        self.BQ = PileupRead.alignment.query_qualities[
            PileupRead.query_position]
        self.MQ = PileupRead.alignment.mapq
        self.is_reverse = PileupRead.alignment.is_reverse


def is_reverse_to_str(boolean):
    if(boolean is True):
        return "reverse"
    elif(boolean is False):
        return "forward"
    else:
        return "unmapped"


class AltAggregateInfo:

    """
    Class which holds summary information for a given alt allele
    at a specific position. Meant to be used as part of the PCInfo
    class as values in the VariantDict.
    All alt alleles in this set of reads should be identical.
    recList must be a list of PRInfo objects.
    """

    def __init__(self, recList, consensus="default",
                 mergedSize="default",
                 totalSize="default",
                 minMQ=0,
                 minBQ=0,
                 contig="default",
                 pos="default",
                 DOC="default",
                 DOCTotal="default"):
        import collections
        if(consensus == "default"):
            raise ThisIsMadness("A consensus nucleotide must be provided.")
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
        # Check that all alt alleles are identical
        self.recList = [rec for rec in recList
                        if rec.MQ >= minMQ and rec.BQ >= minBQ]
        try:
            assert(sum([rec.BaseCall == recList[
                0].BaseCall for rec in recList]) == len(recList))
        except AssertionError:
            # print("recList repr: {}".format(repr(recList)))
            print("Alt alleles: {}".format([i.BaseCall for i in recList]))
            raise ThisIsMadness(
                "AltAggregateInfo requires that all alt alleles agree.")
        self.TotalReads = np.sum([rec.FM for rec in recList])
        self.MergedReads = len(recList)
        self.ReverseMergedReads = np.sum([rec.is_reverse for rec in recList])
        self.ForwardMergedReads = self.MergedReads - self.ReverseMergedReads
        self.ReverseTotalReads = np.sum([
            rec.FM for rec in recList if rec.is_reverse])
        self.ForwardTotalReads = self.TotalReads - self.ReverseTotalReads
        self.AveFamSize = float(self.TotalReads) / self.MergedReads
        self.TotalAlleleDict = {"A": 0, "C": 0, "G": 0, "T": 0}
        self.SumBQScore = sum([rec.BQ for rec in recList])
        self.SumMQScore = sum([rec.MQ for rec in recList])
        try:
            self.AveMQ = float(self.SumMQScore) / len(self.recList)
        except ZeroDivisionError:
            self.AveMQ = 0
        try:
            self.AveBQ = float(self.SumBQScore) / len(self.recList)
        except ZeroDivisionError:
            self.AveBQ = 0
        self.ALT = recList[0].BaseCall
        self.consensus = consensus
        self.minMQ = minMQ
        self.minBQ = minBQ

        # Dealing with transitions (e.g., A-T) and their strandedness
        self.transition = "->".join([consensus, self.ALT])
        self.strandedTransitions = {}
        self.strandedTransitionDict = collections.Counter(
            ["&&".join([self.transition, is_reverse_to_str(
                rec.is_reverse)]) for rec in self.recList])
        self.strandedTotalTransitionDict = collections.Counter(
            ["&&".join([self.transition, is_reverse_to_str(
                rec.is_reverse)]) for rec in self.recList] * rec.FM)
        self.strandedMergedTransitionDict = collections.Counter(
            ["&&".join([self.transition, is_reverse_to_str(
                rec.is_reverse)]) for rec in self.recList])
        self.StrandCountsDict = {}
        self.StrandCountsDict['reverse'] = sum([
            rec.is_reverse for rec in self.recList])
        self.StrandCountsDict['forward'] = sum([
            rec.is_reverse is False for rec in self.recList])
        self.StrandCountsTotalDict = {}
        self.StrandCountsTotalDict['reverse'] = sum(
            [rec.FM for rec in self.recList if rec.is_reverse is True])
        self.StrandCountsTotalDict['forward'] = sum(
            [rec.FM for rec in self.recList if rec.is_reverse is False])
        self.TotalAlleleFrequency = self.TotalReads / float(totalSize)
        self.MergedAlleleFrequency = self.MergedReads / float(mergedSize)
        if(self.ReverseMergedReads != 0 and self.ForwardMergedReads != 0):
            self.BothStrandSupport = True
        else:
            self.BothStrandSupport = False


class PCInfo:

    """
    Takes a pysam.PileupColumn covering one base in the reference
    and makes a new class which has "reference" (inferred from
    consensus) and a list of PRData (one for each read).
    """

    def __init__(self, PileupColumn, minBQ=0, minMQ=0):
        minMQ = int(minMQ)
        minBQ = int(minBQ)
        from collections import Counter
        PysamToChrDict = utilBMF.HTSUtils.GetRefIdDicts()['idtochr']
        self.contig = PysamToChrDict[PileupColumn.reference_id]
        self.pos = PileupColumn.reference_pos
        self.FailedBQReads = sum(
            [pileupRead.alignment.query_qualities
             [pileupRead.query_position] < minBQ
             for pileupRead in PileupColumn.pileups])
        self.FailedMQReads = sum(
            [pileupRead.alignment.mapping_quality < minMQ
             for pileupRead in PileupColumn.pileups])
        self.PCol = PileupColumn
        '''
        if(self.FailedMQReads != 0):
            print("Mapping qualities of reads which failed: {}".format(
                  str([pu.alignment.mapping_quality
                      for pu in PileupColumn.pileups
                      if pu.alignment.mapping_quality >= minMQ])))
        print("Number of reads failed for BQ < {}: {}".format(minBQ,
              self.FailedBQReads))
        print("Number of reads failed for MQ < {}: {}".format(minMQ,
              self.FailedMQReads))
        '''
        self.Records = [
            PRInfo(pileupRead) for pileupRead in PileupColumn.pileups
            if(pileupRead.alignment.mapq >= minMQ) and
            (pileupRead.alignment.query_qualities[
                pileupRead.query_position] >= minBQ)]
        self.MergedReads = len(self.Records)
        try:
            self.TotalReads = sum([rec.FM for rec in self.Records])
        except KeyError:
            self.TotalReads = self.MergedReads
        try:
            self.consensus = Counter(
                [rec.BaseCall for rec in self.Records]).most_common(1)[0][0]
        except IndexError:
            self.consensus = Counter([rec.BaseCall
                                      for rec in
                                      [PRInfo(pileupRead) for pileupRead in
                                       PileupColumn.pileups]]).most_common(
                1)[0][0]
        self.VariantDict = {}
        for alt in list(set([rec.BaseCall for rec in self.Records])):
            self.VariantDict[alt] = [
                rec for rec in self.Records if rec.BaseCall == alt]
        self.AltAlleleData = [AltAggregateInfo(
                              self.VariantDict[key],
                              consensus=self.consensus,
                              mergedSize=self.MergedReads,
                              totalSize=self.TotalReads,
                              minMQ=minMQ,
                              minBQ=minBQ,
                              contig=self.contig,
                              pos=self.pos,
                              DOC=self.MergedReads,
                              DOCTotal=self.TotalReads
                              ) for key in self.VariantDict.keys()]
        self.TotalFracDict = {}
        for alt in self.AltAlleleData:
            self.TotalFracDict[
                alt.ALT] = float(alt.TotalReads) / self.TotalReads
        self.MergedFracDict = {}
        for alt in self.AltAlleleData:
            self.MergedFracDict[
                alt.ALT] = float(alt.MergedReads) / self.MergedReads
        # Generates statistics based on transitions, e.g. ref-->alt
        TransMergedCounts = {}
        for alt in self.AltAlleleData:
            try:
                TransMergedCounts[alt.transition] += alt.MergedReads
            except KeyError:
                TransMergedCounts[alt.transition] = alt.MergedReads
        self.TransMergedCounts = TransMergedCounts
        TransTotalCounts = {}
        for alt in self.AltAlleleData:
            try:
                TransTotalCounts[alt.transition] += alt.TotalReads
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
                self.MergedAlleleFreqDict[key] = self.MergedAlleleDict[
                    key] / float(self.MergedReads)
                self.TotalAlleleFreqDict[key] = self.TotalAlleleDict[
                    key] / float(self.TotalReads)
        except ZeroDivisionError:
            for key in self.MergedAlleleDict:
                self.MergedAlleleFreqDict[key] = 0.
                self.TotalAlleleFreqDict[key] = 0.
        self.MergedAlleleCountStr = "\t".join([
            str(i) for i in [self.MergedAlleleDict['A'],
                             self.MergedAlleleDict['C'],
                             self.MergedAlleleDict['G'],
                             self.MergedAlleleDict['T']]])
        self.TotalAlleleCountStr = "\t".join([
            str(i) for i in [self.TotalAlleleDict['A'],
                             self.TotalAlleleDict['C'],
                             self.TotalAlleleDict['G'],
                             self.TotalAlleleDict['T']]])
        self.MergedAlleleFreqStr = "\t".join([
            str(i) for i in [self.MergedAlleleFreqDict['A'],
                             self.MergedAlleleFreqDict['C'],
                             self.MergedAlleleFreqDict['G'],
                             self.MergedAlleleFreqDict['T']]])
        self.TotalAlleleFreqStr = "\t".join([
            str(i) for i in [self.TotalAlleleFreqDict['A'],
                             self.TotalAlleleFreqDict['C'],
                             self.TotalAlleleFreqDict['G'],
                             self.TotalAlleleFreqDict['T']]])
        # MergedStrandednessRatioDict is the fraction of reverse reads
        # supporting an alternate allele.
        # E.g., if 5 support the alt, but 2 are mapped to the reverse
        # strand, the value is 0.4
        self.MergedStrandednessRatioDict = {"A": 0., "C": 0., "G": 0., "T": 0.}
        self.TotalStrandednessRatioDict = {"A": 0., "C": 0., "G": 0., "T": 0.}
        for alt in self.AltAlleleData:
            self.MergedStrandednessRatioDict[
                alt.ALT] = alt.ReverseMergedReads / float(alt.MergedReads)
            self.TotalStrandednessRatioDict[
                alt.ALT] = alt.ReverseTotalReads / float(alt.TotalReads)
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

    def toString(self, header=False):
        outStr = ""
        if(header is True):
            outStr = ("#Chr\tPos (0-based)\tRef (Consensus)\tAlt\tTotal "
                      "Reads\tMerged Reads\tTotal Allele Frequency\tMerged "
                      "Allele Frequency\tReverse Total Reads\tForward Total"
                      " Reads\tReverse Merged Reads\tForward Merged Reads"
                      "\tFraction Of Total Reads\t"
                      "Fraction Of Merged Reads\tAverage "
                      "Family Size\t"
                      "BQ Sum\tBQ Mean\tMQ Sum\tMQ Mean\n")
        for alt in self.AltAlleleData:
            outStr += '\t'.join([
                str(i) for i in [self.contig,
                                 self.pos,
                                 self.consensus,
                                 alt.ALT,
                                 alt.TotalReads,
                                 alt.MergedReads,
                                 alt.TotalAlleleFrequency,
                                 alt.MergedAlleleFrequency,
                                 alt.StrandCountsTotalDict['reverse'],
                                 alt.StrandCountsTotalDict['forward'],
                                 alt.StrandCountsDict['reverse'],
                                 alt.StrandCountsDict['forward'],
                                 self.TotalFracDict[alt.ALT],
                                 self.MergedFracDict[alt.ALT],
                                 alt.AveFamSize,
                                 alt.SumBQScore,
                                 alt.AveBQ,
                                 alt.SumMQScore,
                                 alt.AveMQ]]) + "\n"
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
        self.str = self.toString()

    def updateWithPileupColumn(self, PileupColumn):
        self.end += 1
        assert self.start < self.end  # Interval must be positive...
        try:
            self.TotalCoverage += sum([int(pu.alignment.opt(
                                      "FM")) for pu in PileupColumn.pileups])
        except KeyError:
            self.TotalCoverage = PileupColumn.nsegments
        self.UniqueCoverage += PileupColumn.nsegments
        self.AvgTotalCoverage = self.TotalCoverage / float(
            self.end - self.start)
        self.AvgUniqueCoverage = self.UniqueCoverage / float(
            self.end - self.start)

    def toString(self):
        PysamToChrDict = HTSUtils.GetRefIdDicts()['idtochr']
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
    1. Counts for all consensus-->alt transitions (By Total and Merged reads)
    2. Counts for the above, specifying strandedness
    3. Number of Merged Reads supporting each allele
    """
    TransTotalDict = {}
    TransMergedDict = {}
    StrandedTransTotalDict = {}
    StrandedTransMergedDict = {}
    NumTransitionsTotal = 0
    NumTransitionsMerged = 0
    MergedReadsProcessed = 0
    TotalReadsProcessed = 0
    if(PileupTsv == "default"):
        PileupTsv = inputBAM[0:-4] + '.Pileup.tsv'
    if(TransitionTable == "default"):
        TransitionTable = inputBAM[0:-4] + '.Trans.tsv'
    if(StrandedTTable == "default"):
        StrandedTTable = inputBAM[0:-4] + '.StrandedTrans.tsv'
    subprocess.check_call(['samtools', 'index', inputBAM])
    if(os.path.isfile(inputBAM + ".bai") is False):
        pl("Couldn't index BAM - coor sorting, then indexing!")
        inputBAM = HTSUtils.CoorSortAndIndexBam(inputBAM, uuid=True)
    NumProcessed = 0  # Number of processed positions in pileup
    bamHandle = pysam.AlignmentFile(inputBAM, "rb")
    PileupHandle = open(PileupTsv, "w")
    TransHandle = open(TransitionTable, "w")
    StrandedTransHandle = open(StrandedTTable, "w")
    FirstLine = True
    pileupIterator = bamHandle.pileup(max_depth=30000)
    while True:
        try:
            pileupColumn = pileupIterator.next()
        except ValueError:
            pl(("Pysam sometimes runs into errors during iteration which"
                " are not handled with any elegance. Continuing!"))
            continue
        except StopIteration:
            pl("Finished iterations.")
            break
        NumProcessed += 1
        if((NumProcessed) % progRepInterval == 0):
            pl("Number of positions processed: {}".format(
                NumProcessed))
            pl("Total reads processed: {}".format(TotalReadsProcessed))
            pl("Merged reads processed: {}".format(MergedReadsProcessed))
        PColSum = PCInfo(pileupColumn, minBQ=minBQ, minMQ=minMQ)
        MergedReadsProcessed += PColSum.MergedReads
        TotalReadsProcessed += PColSum.TotalReads
        if(FirstLine is True):
            PileupHandle.write(PColSum.toString(header=True))
            FirstLine = False
        else:
            PileupHandle.write(PColSum.toString())
        for key in PColSum.TransMergedCounts.keys():
            try:
                TransMergedDict[
                    key] += PColSum.TransMergedCounts[key]
            except KeyError:
                TransMergedDict[
                    key] = PColSum.TransMergedCounts[key]
            NumTransitionsMerged += PColSum.TransMergedCounts[
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
                                                TransTotalDict[key] / float(
                                                    NumTransitionsTotal),
                                                TransMergedDict[key] / float(
                                                    NumTransitionsMerged)))
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
                StrandedTransTotalDict[key] /
                float(NumTransitionsTotal),
                StrandedTransMergedDict[key] /
                float(NumTransitionsMerged),
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
    1. Counts for all consensus-->alt transitions (By Total and Merged reads)
    2. Counts for the above, specifying strandedness
    3. Number of Merged Reads supporting each allele
    """
    TransTotalDict = {}
    TransMergedDict = {}
    StrandedTransTotalDict = {}
    StrandedTransMergedDict = {}
    NumTransitionsTotal = 0
    NumTransitionsMerged = 0
    MergedReadsProcessed = 0
    TotalReadsProcessed = 0
    if(PileupTsv == "default"):
        PileupTsv = inputBAM[0:-4] + '.Pileup.tsv'
    if(TransitionTable == "default"):
        TransitionTable = inputBAM[0:-4] + '.Trans.tsv'
    if(StrandedTTable == "default"):
        StrandedTTable = inputBAM[0:-4] + '.StrandedTrans.tsv'
    if(bedfile == "default"):
        return CustomPileupFullGenome(inputBAM, PileupTsv=PileupTsv,
                                      TransitionTable=TransitionTable,
                                      StrandedTTable=StrandedTTable,
                                      progRepInterval=progRepInterval)
    bedlines = HTSUtils.ParseBed(bedfile)
    try:
        subprocess.check_call(['samtools', 'index', inputBAM])
    except subprocess.CalledProcessError:
        pl("Couldn't index BAM - coor sorting, then indexing!")
        inputBAM = HTSUtils.CoorSortAndIndexBam(inputBAM, uuid=True)
    NumPos = sum([line[2] - line[1] for line in bedlines])
    NumProcessed = 0  # Number of processed positions in pileup
    pl("Number of positions in bed file: {}".format(NumPos))
    bamHandle = pysam.AlignmentFile(inputBAM, "rb")
    PileupHandle = open(PileupTsv, "w")
    TransHandle = open(TransitionTable, "w")
    StrandedTransHandle = open(StrandedTTable, "w")
    FirstLine = True
    for line in bedlines:
        pileupIterator = bamHandle.pileup(line[0], line[1], line[2],
                                          max_depth=30000)
        while True:
            try:
                pileupColumn = pileupIterator.next()
            except ValueError:
                pl(("Pysam sometimes runs into errors during iteration which"
                    " are not handled with any elegance. Continuing!"))
                continue
            except StopIteration:
                pl("Finished iterations.")
                break
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
            if(FirstLine is True):
                PileupHandle.write(PColSum.toString(header=True))
                FirstLine = False
            else:
                PileupHandle.write(PColSum.toString())
            for key in PColSum.TransMergedCounts.keys():
                try:
                    TransMergedDict[
                        key] += PColSum.TransMergedCounts[key]
                except KeyError:
                    TransMergedDict[
                        key] = PColSum.TransMergedCounts[key]
                NumTransitionsMerged += PColSum.TransMergedCounts[
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
                    TransTotalDict[key] /
                    float(NumTransitionsTotal),
                    TransMergedDict[key] /
                    float(NumTransitionsMerged)))
    StrandedTransHandle.write(("Transition+Strandedness\tTotal Reads "
                               "(Unflattened)\tMergedReads With Transition\t"
                               "Fraction Of Total (Unflattened) Transitions"
                               "\tFraction of Merged Transitions\n"))
    for key in StrandedTransTotalDict.keys():
        if(key.split("&&")[0].split("->")[0] != key.split(
                "&&")[0].split("->")[1]):
            StrandedTransHandle.write(
                "{}\t{}\t{}\t{}\t{}\n".format(
                    key,
                    StrandedTransTotalDict[key],
                    StrandedTransMergedDict[key],
                    StrandedTransTotalDict[key] /
                    float(NumTransitionsTotal),
                    StrandedTransMergedDict[key] /
                    float(NumTransitionsMerged),
                ))
    pl("Transition Table: {}".format(TransitionTable))
    pl("Stranded Transition Table: {}".format(StrandedTTable))
    TransHandle.close()
    PileupHandle.close()
    if(CalcAlleleFreq is True):
        AlleleFreqTsv = AlleleFrequenciesByBase(inputBAM)
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
        "AlleleFrequenciesByBase(inputBAM={},".format(inputBAM) +
        " outputTsv={}, progRepInterval=".format(outputTsv) +
        "{})".format(progRepInterval)))
    if(outputTsv == "default"):
        outputTsv = inputBAM[0:-4] + '.allele.freq.tsv'
    try:
        subprocess.check_call(['samtools', 'index', inputBAM])
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
        pileupIterator = inHandle.pileup(max_depth=30000)
        while True:
            try:
                pileup = pileupIterator.next()
            except ValueError:
                pl(("Pysam sometimes runs into errors during iteration which"
                    " are not handled with any elegance. Continuing!"))
                continue
            except StopIteration:
                pl("Finished iterations.")
                break
            NumProcessed += 1
            if(NumProcessed % progRepInterval == 0):
                pl("Number of base positions processed: {}".format(
                    NumProcessed))
            PColInfo = PCInfo(pileup, minMQ=minMQ, minBQ=minBQ)
            outHandle.write(PColInfo.AlleleFreqStr + "\n")
    else:
        bedLines = HTSUtils.ParseBed(bedfile)
        for line in bedLines:
            puIterator = inHandle.pileup(reference=line[0], start=line[1],
                                         end=line[2],
                                         max_depth=30000)
            while True:
                try:
                    pileup = puIterator.next()
                except ValueError:
                    pl(("Pysam sometimes runs into errors during iteration "
                        "which are not handled well. Continuing!"))
                    continue
                except StopIteration:
                    pl("Finished iterations.")
                    break
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


def BamToCoverageBed(inbam, outbed="default", mincov=5, minMQ=0, minBQ=0):
    """
    Takes a bam file and creates a bed file containing each position
    """
    pl(("Command required to reproduce this call: "
        "BamToCoverageBed(\'{}\', outbed=".format(inbam) +
        "\'{}\', mincov={})".format(outbed, mincov)))
    pl(("WARNING: Coverage counts for this script"
        " are wrong. Fix in the works!"
        " It seems to only show up for very long regions."))
    if(outbed == "default"):
        outbed = inbam[0:-4] + ".doc.bed"
    subprocess.check_call(shlex.split("samtools index {}".format(inbam)),
                          shell=False)
    if(os.path.isfile(inbam + ".bai") is False):
        pl("Bam index file was not created. Sorting and indexing.")
        inbam = HTSUtils.CoorSortAndIndexBam(inbam)
        pl("Sorted BAM Location: {}".format(inbam))
    inHandle = pysam.AlignmentFile(inbam, "rb")
    outHandle = open(outbed, "w")
    workingChr = 0
    workingPos = 0
    outHandle.write("\t".join(["#Chr", "Start", "End",
                               "Number Of Merged Mapped Bases",
                               "Number of Unmerged Mapped Bases",
                               "Avg Merged Coverage",
                               "Avg Total Coverage"]) + "\n")
    pl("Beginning PileupToBed.")
    pileupIterator = inHandle.pileup(max_depth=30000)
    ChrToPysamDict = utilBMF.HTSUtils.GetRefIdDicts()['chrtoid']
    while True:
        try:
            PC = PCInfo(pileupIterator.next(), minMQ=minMQ, minBQ=minBQ)
        except ValueError:
            pl(("Pysam sometimes runs into errors during iteration which"
                " are not handled with any elegance. Continuing!"))
            continue
        except StopIteration:
            pl("Finished iterations.")
            break
        if("Interval" not in locals()):
            if(PC.MergedReads >= mincov):
                Interval = PileupInterval(contig=PC.PCol.reference_id,
                                          start=PC.PCol.reference_pos,
                                          end=PC.PCol.reference_pos + 1)
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
                if(PC.PCol.reference_pos - workingPos == 1):
                    try:
                        Interval.updateWithPileupColumn(PC.PCol)
                        workingPos += 1
                    except AssertionError:
                        del Interval
                else:
                    outHandle.write(Interval.toString() + "\n")
                    if(PC.PCol.nsegments >= mincov):
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
                    outHandle.write(Interval.toString() + "\n")
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


def CalcWithinBedCoverage(inbam, bed="default", minMQ=0, minBQ=0,
                          outbed="default"):
    """
    Calculates DOC and creates a bed file containing coverage information
    for each position in a provided bed, only counting reads with
    a given minimum mapping quality or base quality.
    """
    if(outbed == "default"):
        outbed = inbam[0:-4] + ".doc.bed"
    pl(("Command required to reproduce this call: "
        "CalcWithinBedCoverage(\'{}\', bed=".format(inbam) +
        "\'{}\', minMQ=\'{}\', minBQ=".format(bed, minMQ) +
        "\'{}\', outbed={})".format(minBQ, outbed)))
    if(bed == "default"):
        pl("Bed file required for CalcWithinBedCoverage")
        raise ThisIsMadness("Bedfile required for calculating coverage.")
    bedLines = HTSUtils.ParseBed(bed)
    subprocess.check_call(shlex.split("samtools index {}".format(inbam)),
                          shell=False)
    if(os.path.isfile(inbam + ".bai") is False):
        pl("Bam index file could not be created. Sorting and indexing.")
        inbam = HTSUtils.CoorSortAndIndexBam(inbam)
        pl("Sorted BAM Location: {}".format(inbam))
    inHandle = pysam.AlignmentFile(inbam, "rb")
    outHandle = open(outbed, "w")
    outHandle.write("\t".join(["#Chr", "Start", "End",
                               "Number Of Merged Mapped Bases",
                               "Number of Unmerged Mapped Bases",
                               "Avg Merged Coverage",
                               "Avg Total Coverage"]) + "\n")
    for line in bedLines:
        TotalReads = 0
        MergedReads = 0
        pileupIterator = inHandle.pileup(line[0],
                                         line[1],
                                         line[2],
                                         max_depth=30000)
        while True:
            try:
                PC = PCInfo(pileupIterator.next(), minMQ=minMQ, minBQ=minBQ)
            except ValueError:
                pl(("Pysam sometimes runs into errors during iteration which"
                    " are not handled with any elegance. Continuing!"))
                continue
            except StopIteration:
                pl("Finished iterations.")
                break
            TotalReads += PC.TotalReads
            MergedReads += PC.MergedReads
        outHandle.write(
            "\t".join(
                [str(i)
                 for i in
                 [line[0],
                  line[1],
                  line[2],
                  MergedReads, TotalReads, MergedReads /
                  float(line[2] - line[1]),
                  TotalReads / float(line[2] - line[1])]]) + "\n")
    inHandle.close()
    outHandle.close()
    return outbed


def CalcWithoutBedCoverage(inbam, bed="default", minMQ=0, minBQ=0,
                           outbed="default"):
    """
    Calculates DOC and creates a bed file containing each position
    not in provided bed, only counting reads with a given minimum mapping
    quality or base quality.
    """
    if(outbed == "default"):
        outbed = inbam[0:-4] + ".doc.SBI.bed"
    pl(("Command required to reproduce this call: "
        "CalcWithoutBedCoverage(\'{}\', bed=".format(inbam) +
        "\'{}\', minMQ=\'{}\', minBQ=".format(bed, minMQ) +
        "\'{}\', outbed={})".format(minBQ, outbed)))
    if(bed == "default"):
        pl("Bed file required for CalcWithoutBedCoverage")
        raise ThisIsMadness("Bedfile required for calculating coverage.")
    bedLines = HTSUtils.ParseBed(bed)
    subprocess.check_call(shlex.split("samtools index {}".format(inbam)),
                          shell=False)
    if(os.path.isfile(inbam + ".bai") is False):
        pl("Bam index file could not be created. Sorting and indexing.")
        inbam = HTSUtils.CoorSortAndIndexBam(inbam)
        pl("Sorted BAM Location: {}".format(inbam))
    inHandle = pysam.AlignmentFile(inbam, "rb")
    outHandle = open(outbed, "w")
    outHandle.write("\t".join(["#Chr", "Start", "End",
                               "Number Of Merged Mapped Bases",
                               "Number of Unmerged Mapped Bases"]) + "\n")
    pileupIterator = inHandle.pileup(max_depth=30000)
    while True:
        try:
            PC = PCInfo(pileupIterator.next(), minMQ=minMQ, minBQ=minBQ)
            if HTSUtils.PosContainedInBed(PC.contig, PC.pos, bedLines):
                continue
            outHandle.write("\t".join(
                [str(i)
                 for i in
                 [PC.contig, PC.pos, PC.pos + 1, PC.MergedReads,
                  PC.TotalReads]]) + "\n")
        except ValueError:
            pl(("Pysam sometimes runs into errors during iteration which"
                " are not handled with any elegance. Continuing!"))
            continue
        except StopIteration:
            pl("Finished iterations.")
            break
    inHandle.close()
    outHandle.close()
    pass
