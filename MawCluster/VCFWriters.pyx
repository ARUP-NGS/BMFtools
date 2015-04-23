# cython: boundscheck=False, c_string_type=str, c_string_encoding=ascii
# cython: cdivision=True, cdivision_warnings=True, profile=True
import logging
from operator import attrgetter as oag
import sys
from subprocess import CalledProcessError

import pysam
import cython
from cytoolz import map as cmap

from MawCluster.SNVUtils import GetVCFHeader
from MawCluster.PileupUtils import (pPileupColumn,
                                    GetDiscordantReadPairs, PCInfo)
from MawCluster.SNVUtils cimport VCFPos
from utilBMF.HTSUtils import (PysamToChrDict, printlog as pl,
                              ParseBed, PopenDispatcher, PopenCall)
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
               cython.str experiment="",
               cython.bint parallel=False):
    pl("Command to reproduce function call: "
       "SNVCrawler('{}', bed=\"{}\"".format(inBAM, bed) +
       ", minMQ={}, minBQ={}, OutVCF".format(minMQ, minBQ) +
       "=\"{}\", MaxPValue={}".format(OutVCF, MaxPValue) +
       ",keepConsensus={}, reference=".format(keepConsensus) +
       "\"{}\", reference_is_path={}".format(reference, reference_is_path) +
       "commandStr=\"{}\", fileFormat=\"{}\"".format(commandStr, fileFormat) +
       ", FILTERTags=\"{}\", INFOTags=\"{}\"".format(FILTERTags, INFOTags) +
       ", FORMATTags=\"{}\", writeHeader={}".format(FORMATTags, writeHeader) +
       ", minFracAgreed={}, minFA={})".format(minFracAgreed, minFA))
    cdef cython.long NumDiscordantPairs = 0
    cdef cython.str VCFLineString, VCFString
    cdef pysam.calignmentfile.IteratorColumnRegion ICR
    cdef pysam.calignmentfile.IteratorColumnAllRefs ICAR
#   cdef pysam.calignmentfile.AlignmentFile discPairHandle, inHandle
    cdef pysam.cfaidx.FastaFile refHandle
    cdef list line, discReads, VCFLines, bedlines
    cdef pPileupRead_t i
    cdef pysam.calignmentfile.AlignedSegment read

    refHandle = pysam.FastaFile(reference)
    if(bed != "default"):
        pl("Bed file used: {}".format(bed))
        bedSet = True
        bedlines = ParseBed(bed)
    else:
        bedSet = False
    if(isinstance(bed, list)):
        bedlines = bed
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
        pl("Bed file provided - iterating through bed columns")
        if(parallel):
            pl("About to make parallel call. WHOO!")
            dispatcher = GetPopenDispatcherICRs(
                bedlines, inBAM, minMQ=minMQ, minFA=minFA, minBQ=minBQ,
                experiment=experiment, minFracAgreed=minFracAgreed,
                MaxPValue=MaxPValue,
                ref=reference, keepConsensus=keepConsensus, threads=2)
            dispatcher.submit()
            if dispatcher.daemon() == 0:
                ohw("\n".join(dispatcher.outstrs.values()) + "\n")
            else:
                raise CalledProcessError(sum(
                    [d.poll() for d in dispatcher.dispatches],
                    "\n".join([d.commandString for
                               d in dispatcher.dispatches])))
        else:
            for line in bedlines:
                ICR = pileupCall(line[0], line[1],
                                 max_depth=200000,
                                 multiple_iterators=True)
                VCFLines, discReads = IteratorColumnRegionToTuple(
                    ICR, minMQ=minMQ, minBQ=minBQ, minFA=minFA,
                    minFracAgreed=minFracAgreed, experiment=experiment,
                    MaxPValue=MaxPValue, puEnd=line[2], refHandle=refHandle,
                    keepConsensus=keepConsensus, reference=reference)
                ohw("\n".join(VCFLines) + "\n")
    else:
        ICAR = pileupCall(max_depth=200000, multiple_iterators=True)
        PileupIt = ICAR.next
        while True:
            try:
                VCFString, discReads = PileupItToVCFLines(
                    PileupIt(), minMQ=minMQ, minBQ=minBQ, minFA=minFA,
                    minFracAgreed=minFracAgreed, MaxPValue=MaxPValue,
                    experiment=experiment, reference=reference,
                    refHandle=refHandle, keepConsensus=keepConsensus)
                for read in discReads:
                    dpw(read)
            except StopIteration:
                pl("Finished iterations.")
                break
            if(len(VCFString) != 0):
                ohw(VCFString + "\n")
    discPairHandle.close()
    return OutVCF


@cython.returns(tuple)
def PileupItToVCFLines(pysam.calignmentfile.PileupColumn PileupCol,
                       cython.long minMQ=-1, cython.long minBQ=-1,
                       cython.str experiment="",
                       cython.long minFA=-1, cython.float minFracAgreed=-1.,
                       cython.float MaxPValue=-1.,
                       cython.bint keepConsensus=False,
                       cython.str reference="default",
                       pysam.cfaidx.FastaFile refHandle=None):
    cdef pPileupColumn_t PileupColumn
    cdef PCInfo_t PC
    cdef list DiscRPs, reads, discReads, DiscRPNames
    cdef pysam.calignmentfile.AlignedSegment read
    cdef pPileupRead_t i
    cdef VCFPos_t pos
    cdef cython.long NumDiscordantPairs
    if(refHandle is None):
        raise ThisIsMadness("refHandle must be provided to write VCF lines!")
    PileupColumn = pPileupColumn(PileupCol)
    PC = PCInfo(PileupColumn, minMQ=minMQ, minBQ=minBQ, experiment=experiment,
                minFracAgreed=0.0, minFA=2,
                experiment=experiment)
    DiscRPs = GetDiscordantReadPairs(PileupColumn)
    DiscRPNames = list(set(cmap(oag("name"), DiscRPs)))
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
def IteratorColumnRegionToTuple(
        pysam.calignmentfile.IteratorColumnRegion ICR,
        cython.long minMQ=-1, cython.long minBQ=-1, cython.str experiment="",
        cython.long minFA=-1, cython.float minFracAgreed=-1.,
        cython.float MaxPValue=-1., puEnd=-1,
        pysam.cfaidx.FastaFile refHandle=None,
        cython.bint keepConsensus=True,
        cython.str reference="default"):
    cdef pysam.calignmentfile.PileupColumn psPileupColumn
    cdef list VCFLines, discReads, allDiscReads
    if(puEnd < 0):
        raise ThisIsMadness("Need final entry in bed line to do ICR->VCF")
    VCFLines = []
    allDiscReads = []
    for psPileupColumn in ICR:
        posStr, discReads = PileupItToVCFLines(
            psPileupColumn, minMQ=minMQ, minBQ=minBQ, minFA=minFA,
            minFracAgreed=minFracAgreed, experiment=experiment,
            MaxPValue=MaxPValue, refHandle=refHandle,
            keepConsensus=keepConsensus, reference=reference)
        if(len(posStr) != 0):
            VCFLines.append(posStr)
        allDiscReads = list(set(allDiscReads + discReads))
        if(psPileupColumn.pos > puEnd):
            break
    return VCFLines, allDiscReads


@cython.returns(cython.str)
def IteratorColumnRegionToStr(
        pysam.calignmentfile.IteratorColumnRegion ICR,
        cython.long minMQ=-1, cython.long minBQ=-1, cython.str experiment="",
        cython.long minFA=-1, cython.float minFracAgreed=-1.,
        cython.float MaxPValue=-1., puEnd=-1,
        pysam.cfaidx.FastaFile refHandle=None,
        cython.bint keepConsensus=True,
        cython.str reference="default"):
    cdef pysam.calignmentfile.PileupColumn psPileupColumn
    cdef list VCFLines, discReads, allDiscReads
    if(puEnd < 0):
        raise ThisIsMadness("Need final entry in bed line to do ICR->VCF")
    VCFLines = []
    allDiscReads = []
    for psPileupColumn in ICR:
        posStr, discReads = PileupItToVCFLines(
            psPileupColumn, minMQ=minMQ, minBQ=minBQ, minFA=minFA,
            minFracAgreed=minFracAgreed, experiment=experiment,
            MaxPValue=MaxPValue, refHandle=refHandle,
            keepConsensus=keepConsensus, reference=reference)
        if(len(posStr) != 0):
            VCFLines.append(posStr)
        allDiscReads = list(set(allDiscReads + discReads))
        if(psPileupColumn.pos > puEnd):
            break
    return "\n".join(VCFLines)


@cython.returns(pysam.calignmentfile.IteratorColumnRegion)
def GetICR(line, pysam.calignmentfile.AlignmentFile inHandle=None):
    """
    Gets an ICR from an alignment file and a bed line.
    """
    try:
        return inHandle.pileup(line[0], line[1], line[2])
    except IndexError:
        try:
            return inHandle.pileup(line[0], line[1])
        except IndexError:
            raise ThisIsMadness("bed line has only one field?")


def GetICRString(bedline, bampath, minMQ=-1, minFA=-1, minBQ=-1,
                 experiment="", minFracAgreed=-1., MaxPValue=-1.,
                 ref="default", keepConsensus=True):
    if(ref == "default"):
        raise ThisIsMadness("Hey, I need a path to an faidx'd "
                            "fasta reference to variant-call.")
    puEnd = bedline[2]
    string = ("python -c 'import pysam; af = pysam.AlignmentFile("
              "\"%s\");  b = af.pileup" % bampath +
              "(\"%s\", %s, max_depth=200000," % (bedline[0], bedline[1]) +
              " multiple_iterators=True);"
              "from MawCluster.VCFWriters import IteratorColumnRegionToStr;" +
              "c = IteratorColumnRegionToStr(b, minMQ=%s, minBQ=" % (minMQ) +
              "%s, minFA=%s, experiment=\"%s\", minFr" % (minBQ, minFA,
                                                          experiment) +
              "acAgreed=%s, MaxPValue=%s, puEn" % (minFracAgreed, MaxPValue) +
              "d=%s, refHandle=pysam.FastaFile(\"%s\"), " % (puEnd, ref) +
              "keepConsensus=%s, reference=\"%s\"" % (keepConsensus, ref) +
              "); import sys; sys.stdout.write(c)'")
    return string


def GetPopenDispatcherICRs(bed, bampath, minMQ=-1, minFA=-1, minBQ=-1,
                           experiment="", minFracAgreed=-1., MaxPValue=-1.,
                           ref="default", keepConsensus=True, threads=2):
    if isinstance(bed, str):
        bed = ParseBed(bed)
    return PopenDispatcher([GetICRString(bedline, bampath, minMQ=minMQ,
                                         minFA=minFA, minBQ=minBQ,
                                         MaxPValue=MaxPValue, ref=ref,
                                         keepConsensus=keepConsensus,
                                         minFracAgreed=minFracAgreed) for
                            bedline in bed], threads=threads)
