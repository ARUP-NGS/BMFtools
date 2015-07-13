# cython: boundscheck=False, c_string_type=str, c_string_encoding=ascii

# Standard library/builtins.
import logging
from operator import attrgetter as oag
import sys
from subprocess import CalledProcessError, check_call
import warnings

# Third party includes
import pysam
import cython

# local python imports
cimport MawCluster.PileupUtils as PileupUtils
from . import BCVCF
from .SNVUtils import GetVCFHeader
from .PileupUtils import (pPileupColumn,
                          GetDiscordantReadPairs, PCInfo)
from utilBMF.HTSUtils import (printlog as pl,
                              ParseBed, PopenDispatcher, PopenCall,
                              parseConfig, TrimExt, BedListToStr)
from utilBMF import HTSUtils
from utilBMF.ErrorHandling import (ThisIsMadness, FunctionCallException,
                                   AbortMission)
from utilBMF.ArgumentSketcher import parseSketchConfig

# local cython imports
from utilBMF.HTSUtils cimport cystr
from MawCluster.SNVUtils cimport VCFPos
from utilBMF.HTSUtils cimport pPileupRead

# Third party cython imports
from cython cimport str as cystr
cimport cython
cimport pysam.calignmentfile
cimport pysam.cfaidx

# ctypedefs
ctypedef PileupUtils.PCInfo PCInfo_t
ctypedef PileupUtils.pPileupColumn pPileupColumn_t
ctypedef pPileupRead pPileupRead_t
ctypedef VCFPos VCFPos_t


"""
Programs which write VCFs.
Currently: SNVCrawler.

"""


@cython.returns(cystr)
def SNVCrawler(inBAM,
               cystr bed="default",
               int minMQ=0,
               int minBQ=0,
               cystr OutVCF="default",
               float MaxPValue=1e-30,
               cython.bint keepConsensus=False,
               cystr reference="default",
               cython.bint reference_is_path=False,
               cystr commandStr="default",
               cystr fileFormat="default",
               cystr FILTERTags="default",
               cystr INFOTags="default",
               cystr FORMATTags="default",
               cython.bint writeHeader=True,
               float minFracAgreed=0.0,
               int minFA=2,
               cystr experiment="",
               cython.bint parallel=True,
               sampleName="DefaultSampleName",
               conf="default"):
    cdef int NumDiscordantPairs = 0
    cdef cystr VCFLineString, VCFString
    cdef pysam.calignmentfile.IteratorColumnRegion ICR
    cdef pysam.calignmentfile.IteratorColumnAllRefs ICAR
    cdef pysam.calignmentfile.PileupColumn PILEUP_COLUMN
#   cdef pysam.calignmentfile.AlignmentFile discPairHandle, inHandle
    cdef pysam.cfaidx.FastaFile refHandle
    cdef list line, discReads, VCFLines, bedlines
    cdef pPileupRead_t i
    cdef pysam.calignedsegment.AlignedSegment read
    cdef pPileupColumn_t pPC
    confDict = parseSketchConfig(conf) if(conf != "default") else {}
    if(bed != "default"):
        pl("Bed file used: {}".format(bed))
        bedSet = True
        bedlines = ParseBed(bed)
    else:
        if("bed" not in confDict.iterkeys()):
            print("no bed file provided...")
            bedSet = False
        else:
            bedlines = ParseBed(confDict["bed"])
    if(isinstance(bed, list)):
        bedlines = bed
    if("confDict" in locals()):
        try:
            minMQ = int(confDict['minMQ'])
        except KeyError:
            pass
        try:
            minBQ = int(confDict['minBQ'])
        except KeyError:
            pass
        try:
            MaxPValue = float(confDict['MaxPValue'])
        except KeyError:
            pass
        try:
            minFA = int(confDict['minFA'])
        except KeyError:
            pass
        try:
            reference = confDict['ref']
        except KeyError:
            pass
        try:
            minFracAgreed = float(confDict['minFracAgreed'])
        except KeyError:
            pass
    refHandle = pysam.FastaFile(reference)
    if(OutVCF == "default"):
        OutVCF = TrimExt(inBAM) + ".bmf.vcf"
    inHandle = pysam.AlignmentFile(inBAM, "rb")
    if(OutVCF == "stdout"):
        pl("Writing to stdout")
        outHandle = sys.stdout
    else:
        pl("Writing to file %s" % OutVCF)
        outHandle = open(OutVCF, "w")
    pileupCall = inHandle.pileup
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
        pl("Bed file provided - iterating through bed columns")
        for line in bedlines:
            pl("Making pileup call for region %s" % line)
            ICR = pileupCall(line[0], line[1],
                             max_depth=200000,
                             multiple_iterators=False)
            PileupIt = ICR.next
            while True:
                try:
                    PILEUP_COLUMN = PileupIt()
                    while(PILEUP_COLUMN.pos < line[1]):
                        PILEUP_COLUMN = PileupIt()
                    pPC = pPileupColumn(PILEUP_COLUMN)
                except StopIteration:
                    pl("Finishing iterations for bed line: %s" % repr(
                        line))
                    break
                except ValueError:
                    pl("Pysam's heinous errors in iteration happened...")
                    pl("Region: %s" % repr(line))
                    raise ValueError(repr(line))
                posStr = pPileupColToVCFLines(
                    pPC, minMQ=minMQ, minBQ=minBQ, minFA=minFA,
                    minFracAgreed=minFracAgreed, experiment=experiment,
                    MaxPValue=MaxPValue, refHandle=refHandle,
                    keepConsensus=keepConsensus, reference=reference)
                if(posStr != "empty"):
                    ohw(posStr + "\n")
                if(pPC.reference_pos > line[2]):
                    pl("Bed region %s complete." % BedListToStr(line))
                    break
    else:
        ICAR = pileupCall(max_depth=200000, multiple_iterators=False)
        PileupIt = ICAR.next
        while True:
            try:
                VCFString = PileupItToVCFLines(
                    PileupIt(), minMQ=minMQ, minBQ=minBQ, minFA=minFA,
                    minFracAgreed=minFracAgreed, MaxPValue=MaxPValue,
                    experiment=experiment, reference=reference,
                    refHandle=refHandle, keepConsensus=keepConsensus)
            except StopIteration:
                pl("Finished iterations.")
                break
            ohw(VCFString + "\n")
    outHandle.close()
    return OutVCF


@cython.returns(cystr)
def PileupItToVCFLines(pysam.calignmentfile.PileupColumn PileupCol,
                       int minMQ=-1, int minBQ=-1,
                       cystr experiment="",
                       int minFA=-1, float minFracAgreed=-1.,
                       float MaxPValue=-1.,
                       cython.bint keepConsensus=False,
                       cystr reference="default",
                       pysam.cfaidx.FastaFile refHandle=None):
    cdef pPileupColumn_t PileupColumn
    cdef PCInfo_t PC
    cdef list DiscRPs, reads, discReads, DiscRPNames
    cdef pysam.calignedsegment.AlignedSegment read
    cdef pPileupRead_t i
    cdef VCFPos_t pos
    cdef int NumDiscordantPairs
    if(refHandle is None):
        raise ThisIsMadness("refHandle must be provided to write VCF lines!")
    PileupColumn = pPileupColumn(PileupCol)
    PC = PCInfo(PileupColumn, minMQ=minMQ, minBQ=minBQ,
                experiment=experiment,
                minFracAgreed=minFracAgreed, minFA=minFA)
    if(PC.MergedReads == 0):
        print("PC Failure log: %s" % PC.getQCFailString())
        return "empty"
    pos = VCFPos(PC, MaxPValue=MaxPValue,
                 keepConsensus=keepConsensus,
                 reference=reference,
                 minFracAgreed=minFracAgreed,
                 minFA=minFA, refHandle=refHandle,
                 NDP=len(PC.DiscNames))
    if(len(pos.VCFLines) == 0):
        return "empty"
    return str(pos)


@cython.returns(cystr)
def pPileupColToVCFLines(pPileupColumn_t PileupColumn,
                         int minMQ=-1, int minBQ=-1,
                         cystr experiment="",
                         int minFA=-1, float minFracAgreed=-1.,
                         float MaxPValue=-1.,
                         cython.bint keepConsensus=False,
                         cystr reference="default",
                         pysam.cfaidx.FastaFile refHandle=None):
    """
    Returns 'empty' if no positions made it.
    """
    cdef PCInfo_t PC
    cdef VCFPos_t pos
    if(refHandle is None):
        raise ThisIsMadness("refHandle must be provided to write VCF lines!")
    PC = PCInfo(PileupColumn, minMQ=minMQ, minBQ=minBQ,
                experiment=experiment,
                minFracAgreed=minFracAgreed, minFA=minFA,
                experiment=experiment)
    pos = VCFPos(PC, MaxPValue=MaxPValue,
                 keepConsensus=keepConsensus,
                 reference=reference,
                 minFracAgreed=minFracAgreed,
                 minFA=minFA, refHandle=refHandle,
                 NDP=len(PC.DiscNames))
    if(len(pos.VCFLines) == 0):
        return "empty"
    return str(pos)


@cython.returns(cystr)
def IteratorColumnRegionToStr(
        pysam.calignmentfile.IteratorColumnRegion ICR,
        int minMQ=-1, int minBQ=-1, cystr experiment="",
        int minFA=-1, float minFracAgreed=-1.,
        float MaxPValue=-1., puEnd=-1,
        pysam.cfaidx.FastaFile refHandle=None,
        cython.bint keepConsensus=True,
        cystr reference="default"):
    cdef pysam.calignmentfile.PileupColumn psPileupColumn
    cdef list VCFLines, discReads, allDiscReads
    if(puEnd < 0):
        raise ThisIsMadness("Need final entry in bed line to do ICR->VCF")
    VCFLines = []
    allDiscReads = []
    vla = VCFLines.append
    for psPileupColumn in ICR:
        posStr = PileupItToVCFLines(
            psPileupColumn, minMQ=minMQ, minBQ=minBQ, minFA=minFA,
            minFracAgreed=minFracAgreed, experiment=experiment,
            MaxPValue=MaxPValue, refHandle=refHandle,
            keepConsensus=keepConsensus, reference=reference)
        vla(posStr)
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
    cStr = ("python -c 'import pysam; af = pysam.AlignmentFile("
            "\"%s\");  b = af.pileup" % bampath +
            "(\"%s\", %s, max_depth=200000," % (bedline[0], bedline[1]) +
            " multiple_iterators=False);"
            "from MawCluster.VCFWriters import IteratorColumnRegionToStr;" +
            "c = IteratorColumnRegionToStr(b, minMQ=%s, minBQ=" % (minMQ) +
            "%s, minFA=%s, experiment=\"%s\", minFr" % (minBQ, minFA,
                                                        experiment) +
            "acAgreed=%s, MaxPValue=%s, puEn" % (minFracAgreed, MaxPValue) +
            "d=%s, refHandle=pysam.FastaFile(\"%s\"), " % (puEnd, ref) +
            "keepConsensus=%s, reference=\"%s\"" % (keepConsensus, ref) +
            "); import sys; sys.stdout.write(c)'")
    return cStr


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


def PSNVCall(inBAM, conf="default", threads=-1, outVCF="default"):
        config = parseConfig(conf)
        if("threads" in config.keys() and threads < 0):
            threads = int(config["threads"])
        elif(threads < 0):
            threads = 4  # If threads is set in kwargs, keep it.
        if(outVCF == "default"):
            outVCF = TrimExt(inBAM) + ".psnv.bmf.vcf"
        outHandle = open(outVCF, "w")
        outHandle.write(GetVCFHeader(reference=config["ref"],
                        header=pysam.AlignmentFile(inBAM, "rb").header))
        pl("Splitting BAM file by contig.")
        Dispatcher = HTSUtils.GetBMFsnvPopen(inBAM, config['bed'],
                                             conf=conf,
                                             threads=threads)
        if(Dispatcher.daemon() != 0):
            raise ThisIsMadness("Dispatcher failed somehow.")
        pl("Shell calls completed without errors.")
        for commandString, vcffile in Dispatcher.outstrs.items():
            if(vcffile is None):
                raise FunctionCallException(
                    commandString,
                    "Attempt to cat this vcf to final file failed.", -1)
            check_call("cat %s | grep -v '^#' >> %s" % (vcffile, outVCF),
                       shell=True)
        pl("Filtering VCF by bed file. Pre-filter path: %s" % outVCF)
        bedFilteredVCF = BCVCF.FilterVCFFileByBed(
                    outVCF, bedfile=config["bed"])
        pl("Filtered VCF: %s" % bedFilteredVCF)
        return bedFilteredVCF
