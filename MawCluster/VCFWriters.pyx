# cython: boundscheck=False, c_string_type=str, c_string_encoding=ascii
# cython: cdivision=True, cdivision_warnings=True, profile=True
import logging
from operator import attrgetter as oag
import sys

import pysam
import cython
cimport cython
from cytoolz import map as cmap
cimport pysam.calignmentfile
cimport pysam.cfaidx

from MawCluster.SNVUtils import *
from MawCluster.PileupUtils import pPileupColumn, GetDiscordantReadPairs
cimport MawCluster.PileupUtils as PileupUtils
from MawCluster.SNVUtils cimport VCFPos
from utilBMF.HTSUtils import PysamToChrDict
from utilBMF.HTSUtils cimport pPileupRead
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


@cython.locals(MaxPValue=cython.float, reference_is_path=cython.bint,
               writeHeader=cython.bint, minFA=cython.long)
def SNVCrawler(inBAM,
               bed="default",
               minMQ=0,
               minBQ=0,
               OutVCF="default",
               MaxPValue=1e-30,
               keepConsensus=False,
               reference="default",
               reference_is_path=False,
               commandStr="default",
               fileFormat="default",
               FILTERTags="default",
               INFOTags="default",
               FORMATTags="default",
               writeHeader=True,
               minFracAgreed=0.0, minFA=2,
               experiment=""):
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
    cdef PCInfo_t PC
    cdef pPileupColumn_t PileupColumn
    cdef pysam.calignmentfile.IteratorColumnRegion puIterator
    cdef pysam.calignmentfile.AlignmentFile discPairHandle
    cdef pysam.cfaidx.FastaFile refHandle
    cdef VCFPos_t pos
    cdef list line, DiscRPs, reads
    cdef pPileupRead_t i, read

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
            pl("Combing through bed region {}".format(line),
               level=logging.DEBUG)
            puIterator = pileupCall(line[0], line[1],
                                    max_depth=200000,
                                    multiple_iterators=True)
            PileupIt = puIterator.next
            while True:
                try:
                    PileupColumn = pPileupColumn(PileupIt())
                except StopIteration:
                    pl("Finished iterations.")
                    break
                DiscRPs = GetDiscordantReadPairs(PileupColumn)
                DiscRPNames = list(set(cmap(oag("name"), DiscRPs)))
                PileupColumn.pileups = [i for i in PileupColumn.pileups if
                                        i.alignment.query_name not in DiscRPNames]
                for RP in DiscRPs:
                    reads = RP.RP.getReads()
                    for read in reads:
                        read.set_tag("DP", RP.discordanceString, "Z")
                        discPairHandle.write(read)
                NumDiscordantPairs = len(DiscRPNames)
                PC = PCInfo(PileupColumn, minMQ=minMQ, minBQ=minBQ,
                            experiment=experiment)
                #  pl("Position for pileup (0-based): {}".format(PC.pos),
                #     level=logging.DEBUG)
                if(line[2] <= PC.pos):
                    break
                pos = VCFPos(PC, MaxPValue=MaxPValue,
                             keepConsensus=keepConsensus,
                             reference=reference, refHandle=refHandle,
                             minFracAgreed=minFracAgreed,
                             minFA=minFA, NDP=NumDiscordantPairs)
                VCFLineString = str(pos)
                if(len(VCFLineString) != 0):
                    ohw(VCFLineString + "\n")
                else:
                    logStr = ("Failed to write line!" + repr(dir(pos)) +
                              ("bed line: %s" % line) +
                              "\tref id: %s" % PileupColumn.reference_id +
                              "\tref pos: %s" % PileupColumn.reference_pos)
                    logStr += "AC: %s" % PC.MergedAlleleFreqDict
                    logStr += "PC Str: %s" % str(PC)
                    pl(logStr, level=logging.DEBUG)
                    pl("VCF line not written at position" +
                       PysamToChrDict[PileupColumn.reference_id] + ":" +
                       str(PileupColumn.reference_pos + 1) +
                       " - usually because all reads failed filters.",
                       level=logging.INFO)
    else:
        puIterator = pileupCall(max_depth=200000, multiple_iterators=True)
        PileupIt = puIterator.next
        while True:
            try:
                # Last command - 0 means iterator was where it crashed.
                PileupColumn = pPileupColumn(PileupIt())
                # Last command - 0 means the PCInfo call was where it crashed
                PC = PCInfo(PileupColumn, minMQ=minMQ, minBQ=minBQ,
                            experiment=experiment)
            except ValueError:
                raise ThisIsMadness("Trying to figure out what's going on in "
                                    "pysam's iterations.")
            except StopIteration:
                break
            # TODO: Check to see if it speeds up to not assign and only write.
            DiscRPs = GetDiscordantReadPairs(PileupColumn)
            DiscRPNames = list(set(cmap(oag("name"), DiscRPs)))
            PileupColumn.pileups = [i for i in PileupColumn.pileups if
                                    i.alignment.query_name not in DiscRPNames]
            for RP in DiscRPs:
                reads = RP.RP.getReads()
                for read in reads:
                    read.set_tag("DP", RP.discordanceString, "Z")
                    discPairHandle.write(read)
            NumDiscordantPairs = len(DiscRPNames)
            pos = VCFPos(PC, MaxPValue=MaxPValue,
                         keepConsensus=keepConsensus,
                         reference=reference,
                         minFracAgreed=minFracAgreed,
                         minFA=minFA, refHandle=refHandle,
                         NDP=NumDiscordantPairs)
            VCFLineString = str(pos)
            if(len(VCFLineString) != 0):
                ohw(VCFLineString + "\n")
    discPairHandle.close()
    return OutVCF


# Trying to "parallelize" this...
# I'll get around to it later.
"""


def SNVMinion(inBAM,
              bed="default",
              minMQ=0,
              minBQ=0,
              VCFLines="default",
              MaxPValue=1e-30,
              keepConsensus=False,
              reference="default",
              reference_is_path=False,
              commandStr="default",
              fileFormat="default",
              FILTERTags="default",
              INFOTags="default",
              FORMATTags="default"):
    if(isinstance(bed, str) and bed != "default"):
        bed = HTSUtils.ParseBed(bed)
    if(VCFLines == "default"):
        VCFLines = inBAM[0:-4] + ".bmf.vcf"
    inHandle = pysam.AlignmentFile(inBAM, "rb")
    outHandle = open(VCFLines, "w+")
    if(bed != "default"):
        for line in bed:
            puIterator = pileupCall(line[0], line[1],
                                         max_depth=30000,
                                         multiple_iterators=True)
            while True:
                try:
                    PileupColumn = puIterator.next()
                    PC = PCInfo(PileupColumn, minMQ=minMQ, minBQ=minBQ)
                    # print(PC.__str__())
                except ValueError:
                    pl(("Pysam sometimes runs into errors during iteration w"
                        "hich are not handled with any elegance. Continuing!"))
                    continue
                except StopIteration:
                    pl("Finished iterations.")
                    break
                if(line[2] <= PC.pos):
                    break
                VCFLineString = VCFPos(PC, MaxPValue=MaxPValue,
                                       keepConsensus=keepConsensus,
                                       reference=reference
                                       ).__str__()
                if(len(VCFLineString) != 0):
                    ohw(VCFLineString + "\n")
    else:
        puIterator = pileupCall(max_depth=30000)
        while True:
            try:
                PC = PCInfo(puIterator.next(), minMQ=minMQ, minBQ=minBQ)
            except ValueError:
                pl(("Pysam sometimes runs into errors during iteration which"
                    " are not handled with any elegance. Continuing!"))
                continue
            except StopIteration:
                break
            # TODO: Check to see if it speeds up to not assign and only write.
            VCFLineString = VCFPos(PC, MaxPValue=MaxPValue,
                                   keepConsensus=keepConsensus,
                                   reference=reference).__str__()
            if(len(VCFLineString) != 0):
                ohw(VCFLineString + "\n")
    return VCFLines


def CallSNVCrawler():
    pass


def SNVMaster(inBAM,
              bed="default",
              minMQ=0,
              minBQ=0,
              VCFLines="default",
              MaxPValue=1e-30,
              keepConsensus=False,
              reference="default",
              reference_is_path=False,
              commandStr="default",
              fileFormat="default",
              FILTERTags="default",
              INFOTags="default",
              FORMATTags="default",
              ByContig=True,
              children=2):
    from subprocess import Popen
    if(isinstance(bed, str) and bed != "default"):
        bed = HTSUtils.ParseBed(bed)
    if(VCFLines == "default"):
        VCFLines = inBAM[0:-4] + ".bmf.vcf"
    inHandle = pysam.AlignmentFile(inBAM, "rb")
    outHandle = open(VCFLines, "w")
    ohw(GetVCFHeader(fileFormat=fileFormat,
                                 FILTERTags=FILTERTags,
                                 commandStr=commandStr,
                                 reference=reference,
                                 reference_is_path=False,
                                 header=inHandle.header,
                                 INFOTags=INFOTags,
                                 FORMATTags=FORMATTags))
    if(ByContig):
        contigList = list(set([line[0] for line in bed]))
        jobList = []
        for thread in range(int(children)):
            jobList.append(CallSNVCrawler(inBAM,
                           bed="default",
                           minMQ=0,
                           minBQ=0,
                           VCFLines="default",
                           MaxPValue=1e-30,
                           keepConsensus=False,
                           reference="default",
                           reference_is_path=False,
                           commandStr="default",
                           fileFormat="default",
                           FILTERTags="default",
                           INFOTags="default",
                           FORMATTags="default"))
    pass
"""
