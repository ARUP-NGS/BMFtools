# cython: boundscheck=False, c_string_type=str, c_string_encoding=ascii
# cython: cdivision=True, cdivision_warnings=True, profile=True
import logging
from operator import attrgetter as oag
import sys

import pysam
import cython
from cytoolz import map as cmap

from MawCluster.SNVUtils import GetVCFHeader
from MawCluster.PileupUtils import (pPileupColumn,
                                    GetDiscordantReadPairs, PCInfo)
from MawCluster.SNVUtils cimport VCFPos
from utilBMF.HTSUtils import PysamToChrDict, printlog as pl
from utilBMF.HTSUtils cimport pPileupRead
from utilBMF import HTSUtils
from utilBMF.ErrorHandling import ThisIsMadness
cimport cython
cimport pysam.calignmentfile
cimport pysam.cfaidx
cimport MawCluster.PileupUtils as PileupUtils
ctypedef PileupUtils.PCInfo PCInfo_t
ctypedef PileupUtils.pPileupColumn pPileupColumn_t
ctypedef pPileupRead pPileupRead_t
ctypedef VCFPos VCFPos_t


"""
Programs which write VCFs.
Currently: SNVCrawler.

In development: SV

Settled on two major filters for inclusion in a pileup
1. Min FA (number of family members agreeing) [int=2] (I just can't get
   too big a family...)
2. Min Fraction Agreement Within Family [float=0.6667]

More thoughts?
"""


@cython.returns(cython.str)
def SNVCrawler(inBAM,
               cython.str bed="default",
               cython.long minMQ=0,
               cython.long minBQ=0,
               cython.str OutVCF="default",
               cython.float MaxPValue=1e-30,
               cython.bint keepConsensus=False,
               cython.str reference="default",
               cython.bint reference_is_path=False,
               cython.str commandStr="default",
               cython.str fileFormat="default",
               cython.str FILTERTags="default",
               cython.str INFOTags="default",
               cython.str FORMATTags="default",
               cython.bint writeHeader=True,
               cython.float minFracAgreed=0.0, cython.long minFA=2,
               cython.str experiment=""):
    pl("Command to reproduce function call: "
       "SNVCrawler({}, bed=\"{}\"".format(inBAM, bed) +
       ", minMQ={}, minBQ={}, OutVCF".format(minMQ, minBQ) +
       "=\"{}\", MaxPValue={}".format(OutVCF, MaxPValue) +
       ",keepConsensus={}, reference=".format(keepConsensus) +
       "\"{}\", reference_is_path={}".format(reference, reference_is_path) +
       "commandStr=\"{}\", fileFormat=\"{}\"".format(commandStr, fileFormat) +
       ", FILTERTags=\"{}\", INFOTags=\"{}\"".format(FILTERTags, INFOTags) +
       ", FORMATTags=\"{}\", writeHeader={}".format(FORMATTags, writeHeader) +
       ", minFracAgreed={}, minFA={})".format(minFracAgreed, minFA))
    cdef cython.long NumDiscordantPairs = 0
    cdef cython.str VCFLineString
    cdef pysam.calignmentfile.IteratorColumnRegion ICR
    cdef pysam.calignmentfile.IteratorColumnAllRefs ICAR 
    cdef pysam.calignmentfile.AlignmentFile discPairHandle, inHandle
    cdef pysam.cfaidx.FastaFile refHandle
    cdef list line, discReads, VCFLines
    cdef pPileupRead_t i
    cdef pysam.calignmentfile.AlignedSegment read

    refHandle = pysam.FastaFile(reference)
    if(bed != "default"):
        pl("Bed file used: {}".format(bed))
        bedSet = True
    else:
        bedSet = False
    if(isinstance(bed, str) and bed != "default"):
        bed = HTSUtils.ParseBed(bed)
    if(OutVCF == "default"):
        OutVCF = inBAM[0:-4] + ".bmf.vcf"
    inHandle = pysam.AlignmentFile(inBAM, "rb")
    if(OutVCF == "stdout"):
        outHandle = sys.stdout
    else:
        outHandle = open(OutVCF, "w")
    pileupCall = inHandle.pileup
    discPairHandle = pysam.AlignmentFile(
        inBAM[0:-4] + ".discReadPairs.bam", "wb", template=inHandle)
    ohw = outHandle.write
    dpw = discPairHandle.write
    if(writeHeader):
        try:
            ohw(GetVCFHeader(fileFormat=fileFormat, FILTERTags=FILTERTags,
                             commandStr=commandStr, reference=reference,
                             reference_is_path=False, header=inHandle.header,
                             INFOTags=INFOTags, FORMATTags=FORMATTags))
        except ValueError:
            pl("Looks like the RG header wasn't parseable by pysam - that's u"
               "sually an artefact of the clash between pysam and GATK's ways"
               "of working with RG fields.", level=logging.DEBUG)
            ohw(GetVCFHeader(fileFormat=fileFormat,
                             FILTERTags=FILTERTags,
                             commandStr=commandStr,
                             reference=reference,
                             reference_is_path=False,
                             INFOTags=INFOTags,
                             FORMATTags=FORMATTags))
    if(bedSet):
        for line in bed:
            ICR = pileupCall(line[0], line[1],
                             max_depth=200000,
                             multiple_iterators=True)
            VCFLines, discReads = IteratorColumnRegionToVCFLines(
                ICR, minMQ=minMQ, minBQ=minBQ, minFA=minFA,
                minFracAgreed=minFracAgreed, experiment=experiment,
                MaxPValue=MaxPValue, puEnd=line[2])
            for read in discReads:
                dpw(read)
            ohw("\n".join(VCFLines) + "\n")
    else:
        ICAR = pileupCall(max_depth=200000, multiple_iterators=True)
        PileupIt = puIterator.next
        while True:
            try:
                VCFLineString, discReads = PileupItToVCFLines(
                    PileupIt(), minMQ=minMQ, minBQ=minBQ, minFA=minFA, minFracAgreed=minFracAgreed,
                    MaxPValue=MaxPValue, experiment=experiment)
                for read in discReads:
                    dpw(read)
            except StopIteration:
                    pl("Finished iterations.")
                break
            if(len(VCFLineString) != 0):
                ohw(VCFLineString + "\n")
    discPairHandle.close()
    return OutVCF


@cython.returns(tuple)
def PileupItToVCFLines(pysam.calignmentfile.PileupColumn PileupCol,
                       cython.long minMQ=-1, cython.long minBQ=-1,
                       cython.str experiment="",
                       cython.long minFA=-1, cython.float minFracAgreed=-1.,
                       cython.float MaxPValue=-1.):
    cdef pPileupColumn_t PileupColumn
    cdef PCInfo_t PC
    cdef list DiscRPs, reads, discReads
    cdef pysam.calignmentfile.AlignedSegment read
    cdef pPileupRead_t i
    cdef VCFPos_t pos
    PileupColumn = pPileupColumn(PileupCol)
    PC = PCInfo(PileupColumn, minMQ=minMQ, minBQ=minBQ, experiment=experiment,
                minFracAgreed=0.0, minFA=2,
                experiment=experiment)
    DiscRPs = GetDiscordantReadPairs(PileupColumn)
    PileupColumn.pileups = [i for i in PileupColumn.pileups if
                            i.alignment.query_name not in DiscRPNames]
    discReads = []
    for RP in DiscRPs:
        reads = RP.RP.getReads()
        for read in reads:
            read.set_tag("DP", RP.discordanceString, "Z")
            discReads.append(read)
    NumDiscordantPairs = len(DiscRPNames)
    pos = VCFPos(PC, MaxPValue=MaxPValue,
                 keepConsensus=keepConsensus,
                 reference=reference,
                 minFracAgreed=minFracAgreed,
                 minFA=minFA, refHandle=refHandle,
                 NDP=NumDiscordantPairs)
    return str(pos), discReads


@cython.returns(tuple)
def IteratorColumnRegionToVCFLines(pysam.calignmentfile.IteratorColumnRegion ICR,
                                   cython.long minMQ=-1, cython.long minBQ=-1,
                                   cython.str experiment="",
                                   cython.long minFA=-1, cython.float minFracAgreed=-1.,
                                   cython.float MaxPValue=-1., puEnd=-1):
    cdef pysam.calignmentfile.PileupColumn psPileupColumn
    cdef list VCFLines, discReads, allDiscReads
    if(puEnd < 0):
        raise ThisIsMadness("Need final entry in bed line to do ICR->VCF")
    VCFLines = []
    allDiscReads = []
    for psPileupColumn in ICR:
        posStr, discReads = PileupItToVCFLines(psPileupColumn, minMQ=minMQ, minBQ=minBQ,
                                               minFA=minFA, minFracAgreed=minFracAgreed,
                                               experiment=experiment, MaxPValue=MaxPValue)
        if(len(posStr) != 0)
            VCFLines.append(posStr):
        allDiscReads = list(set(allDiscReads + discReads))
        if(psPileupColumn.pos < puEnd):
            break
    return VCFLines, allDiscReads