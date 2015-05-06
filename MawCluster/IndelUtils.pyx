# cython: c_string_type=str, c_string_encoding=ascii
# cython: profile=True, cdivision=True, cdivision_warnings=True
import operator
import pysam
import cython
from cytoolz import map as cmap
from utilBMF.HTSUtils import (FractionAligned, FractionSoftClipped,
                              printlog as pl, CoorSortAndIndexBam,
                              GetKmersToCheck, FastqStrFromKmerList,
                              BowtieFqToStr, GetMQPassReads,
                              GetInsertionFromAlignedSegment,
                              GetDeletionFromAlignedSegment, IndelQuiver,
                              FacePalm)
from .SNVUtils import (HeaderInfoLine, HeaderFormatLine,
                       HeaderContigLine, HeaderCommandLine,
                       HeaderReferenceLine, HeaderFileFormatLine,
                       GetContigHeaderLines, HeaderFilterLine)
from utilBMF.ErrorHandling import ThisIsMadness
import os.path
import uuid
from cytoolz.itertoolz import frequencies as cyfreq
import sys
cimport pysam.calignmentfile
cimport cython
cimport numpy as np
cimport utilBMF.HTSUtils
ctypedef utilBMF.HTSUtils.IndelQuiver IndelQuiver_t
ctypedef utilBMF.HTSUtils.AbstractIndelContainer AbstractIndelContainer_t
ctypedef np.longdouble_t dtype128_t

"""
Contains utilities for working with indels for HTS data.
"""


def FilterByIndelRelevance(inBAM, indelOutputBAM="default",
                           otherOutputBAM="default",
                           minFamSize=2):
    """
    Writes reads potentially relevant to an indel to indelOutputBAM and
    other reads to the otherOutputBAM.
    idRel stands for indel relevant.
    idIrl stands for indel irrelevant.
    Input BAM must be name sorted: coordinate sorted is not supported!
    """
    cdef cython.long idrel
    cdef cython.long idirl
    cdef pysam.calignmentfile.AlignedSegment read1, read2, entry
    if(indelOutputBAM == "default"):
        indelOutputBAM = ".".join(inBAM.split(".")[0:-1] + ["indRel", "bam"])
    if(otherOutputBAM == "default"):
        otherOutputBAM = ".".join(inBAM.split(".")[0:-1] + ["indIrl", "bam"])
    inHandle = pysam.AlignmentFile(inBAM, "rb")
    indelHandle = pysam.AlignmentFile(indelOutputBAM, "wb", template=inHandle)
    otherHandle = pysam.AlignmentFile(otherOutputBAM, "wb", template=inHandle)
    ohw = otherHandle.write
    ihw = indelHandle.write
    idrel = 0
    idirl = 0
    for entry in inHandle:
        if entry.is_read1:
            read1 = entry
            continue
        else:
            read2 = entry
        assert read1.query_name == read2.query_name
        if(IsIndelRelevant(read1, minFam=minFamSize) or
           IsIndelRelevant(read2, minFam=minFamSize)):
            ihw(read1)
            ihw(read2)
            idrel += 1
        else:
            ohw(read1)
            ohw(read2)
            idirl += 1
    inHandle.close()
    otherHandle.close()
    indelHandle.close()
    pl("Finished filtering by indel relevance. Relevant pairs: %s" % idrel +
       ". Irrelevant pairs: %s" % idirl)
    return indelOutputBAM, otherOutputBAM


def IsIndelRelevant(
        pysam.calignmentfile.AlignedSegment read, cython.long minFam=2,
        cython.float minSF=0.5, cython.bint keepUnmapped=False):
    """
    True if considered relevant to indels.
    False otherwise.
    """
    if(read.opt("FM") < minFam):
        return False
    if(read.cigarstring is None):
        # This read is simply unmapped. Let's give it a chance!
        if(keepUnmapped):
            return True
        return False
    if("I" in read.cigarstring or "D" in read.cigarstring):
        return True
    if(FractionSoftClipped(read) >= minSF):
        return True
    return False


def GetSCFractionArray(inBAM):
    cdef pysam.calignmentfile.AlignedSegment i
    inHandle = pysam.AlignmentFile(inBAM, "rb")
    return [FractionSoftClipped(i) for i in inHandle]


def GetFreebayesCallStr(inBAM, ref="default", bed="default",
                        outVCF="default", ploidy=-1,
                        minMQ=1, minBQ=3, haplotypeLength=80,
                        rdf=1.0, bestNAlleles=100):
    """
    Used to call freebayes for indel calling.
    """
    if(os.path.isfile(inBAM + ".bai") is False):
        pl("BAM must be coordinate-sorted and "
           "indexed for freebayes to call variants. Doing so!")
        inBAM = CoorSortAndIndexBam(inBAM)
    if(outVCF == "default"):
        outVCF = ".".join(inBAM.split(".")[:-1]) + ".fb.vcf"
    if(ref == "default"):
        raise ThisIsMadness("Reference required for freebayes call.")
    if(bed == "default"):
        raise ThisIsMadness("Bed file require for freebayes call currently. "
                            "If there is sufficient demand for an all-regions"
                            " call, I can add that option.")
    if(ploidy < 0):
        pl("Ploidy less than 0. Defaulting to report best N alleles.")
        cStr = ("freebayes -v %s -t %s -f %s " % (outVCF, bed, ref) +
                "-q %s -m %s --haplotype-length" % (minBQ, minMQ) +
                " %s -D %s " % (haplotypeLength, rdf) +
                "--min-alternate-fraction 0 --pooled-continuous" +
                " %s" % inBAM)
    else:
        pl("Ploidy set: %s" % ploidy)
        cStr = ("freebayes -v %s -t %s -f %s " % (outVCF, bed, ref) +
                "-q %s -m %s --haplotype-length" % (minBQ, minMQ) +
                " %s -D %s " % (haplotypeLength, rdf) +
                "--min-alternate-fraction 0 --pooled-continuous" +
                " -p %s --use-best-n-alleles %s" % (ploidy, bestNAlleles) +
                " %s" % inBAM)
    return cStr


@cython.returns(cython.str)
def GetFBOutVCFFromStr(cython.str cStr):
    """
    Gets out vcf from freebayes call. Used for parallelization.
    """
    return cStr.split(" ")[2]


@cython.returns(list)
def GetUniquelyMappableKmers(cython.str ref, cython.long k=30,
                             list bedline=[], cython.long minMQ=1,
                             cython.long padding=-1, cython.long mismatches=2):
    """
    Uses a set of HTSUtils methods to find kmers from a region
    which are uniquely mappable. This makes it possible to do alignment-free
    variant-calling. (Well, except for the bwasw step).
    If no outfile is specified, defaults to stdout
    """
    cdef list kmerList, PassingReadNames
    cdef cython.str fqStr, bowtieStr
    pl("Getting potential kmers")
    kmerList = GetKmersToCheck(ref, k=k, bedline=bedline, padding=padding)
    pl("Making dummy fastq records for each kmer")
    fqStr = FastqStrFromKmerList(kmerList)
    pl("Aligning these kmers to the genome to test for unique mappability"
       " with a given number of mismatches %s and minMQ %s." % (mismatches,
                                                                minMQ))
    bowtieStr = BowtieFqToStr(fqStr, ref=ref, mismatches=mismatches, seed=k)
    PassingReadNames = GetMQPassReads(bowtieStr, minMQ=minMQ)
    return PassingReadNames


@cython.returns(IndelQuiver_t)
def GetTaggedIndelsQuiver(inBAM, cython.long minPairs=2,
                          cython.float minShen=0.2,
                          cython.long window=16, cython.str ref=None,
                          list bedRegion=[], cython.long minMQ=1):
    """
    To PileupTaggedIndels, the bam must have the DSI and DSD tags set
    for reads, and it is strongly recommended that you parse out
    only those reads for this tool.
    """
    Quiver = IndelQuiver(ref=ref, window=window)
    cdef pysam.calignmentfile.AlignedSegment rec
    refHandle = pysam.FastaFile(ref)
    inHandle = pysam.AlignmentFile(inBAM, "rb")
    inFetch = inHandle.fetch(bedRegion[0], bedRegion[1])
    for rec in inFetch:
        if(rec.mapping_quality < minMQ):
            continue
            """
            Since we're looking for DS[ID], if one read has a mapq
            sufficient for inclusion, we get the benefit of seeing that indel
            even if one of the reads failed the MQ threshold.
            """
        try:
            svtags = rec.opt("SV")
        except KeyError:
            FacePalm("SV tags must be marked to call tagged indels!")
        if("DSI" in svtags):
            Quiver.addIndel(GetInsertionFromAlignedSegment(rec, handle=refHandle))
        if("DSD" in svtags):
            Quiver.addIndel(GetDeletionFromAlignedSegment(rec, handle=refHandle))
        if(rec.reference_pos >= bedRegion[2]):
            pl("Finished bed region %s." % bedRegion)
            break
    return Quiver

class IDVCFLine:

    """
    Contains data for writing to a VCF Line, where each ALT has its own line.

    minNumSS is the minumum number of start/stop combinations required to
    support a variant call.
    """
    def __init__(self,
                 AbstractIndelContainer_t IC,
                 cython.str ID=".",
                 cython.long minNumSS=2,
                 cython.str REF=None,
                 dtype128_t reverseStrandFraction=-1.0,
                 cython.bint requireDuplex=True,
                 cython.long minDuplexPairs=2,
                 dtype128_t minFracAgreedForFilter=0.666,
                 cython.long minFA=0, cython.long BothStrandAlignment=-1,
                 dtype128_t pValBinom=0.05,
                 cython.float minAF=-1.,
                 cython.long maxND=10):
        if(BothStrandAlignment < 0):
            raise ThisIsMadness("BothStrandAlignment required for SNVCFLine,"
                                " as it is used in determining whether or no"
                                "t to remove a call for a variant mapping to"
                                " only one strand.")
        self.REF = REF
        self.CHROM = IC.contig
        self.NumStartStops = len(set(IC.StartStops))
        """
        self.NumStartStops = len(set(map(operator.attrgetter("ssString"),
                                         AlleleAggregateObject.recList)))
        This object is in progress.

        self.POS = AlleleAggregateObject.pos + 1
        self.CONS = AlleleAggregateObject.consensus
        self.ALT = AlleleAggregateObject.ALT
        self.QUAL = 1. * AlleleAggregateObject.SumBQScore
        self.BothStrandSupport = AlleleAggregateObject.BothStrandSupport
        self.reverseStrandFraction = reverseStrandFraction
        self.AABothStrandAlignment = BothStrandAlignment
        self.ID = ID
        AC = AlleleAggregateObject.MergedReads
        DOC = AlleleAggregateObject.DOC
        minAAF, maxAAF = ConfidenceIntervalAAF(AC, DOCMerged, pVal=pValBinom)
        try:
            if(float(MaxPValue) < 10 ** (self.QUAL / -10.)):
                if(self.FILTER is not None):
                    self.FILTER += ";LowQual"
                else:
                    self.FILTER = "LowQual"
                self.QUAL *= 0.1
                # This is entirely arbitrary...
        except TypeError:
            print("TypeError! MaxPValue is: {}".format(MaxPValue))
            raise TypeError("MaxPValue not properly parsed.")
        if(AlleleAggregateObject.BothStrandSupport is False and
           BothStrandAlignment):
            if(self.FILTER is not None):
                self.FILTER += ";OneStrandSupport"
            else:
                self.FILTER = "OneStrandSupport"
        if(self.FILTER is None):
            self.FILTER = "PASS"
        if(self.ALT == AlleleAggregateObject.consensus):
            if(self.FILTER is None):
                self.FILTER = "CONSENSUS"
            else:
                self.FILTER += ";CONSENSUS"
        shenRef = min([shen(flankingBefore + REF), shen(REF + flankingAfter)])
        shenVar = min([shen(flankingBefore + self.ALT),
                       shen(self.ALT + flankingAfter)])
        self.InfoFields = {"AC": AC,
                           "AF": 1. * AC / DOC,
                           "BNP": int(-10 * mlog10(pValBinom)),
                           "BSS": AlleleAggregateObject.BothStrandSupport,
                           "TF": 1. * AlleleAggregateObject.TotalReads /
                           AlleleAggregateObject.DOCTotal,
                           "NSS": self.NumStartStops,
                           "MBP": AlleleAggregateObject.MBP,
                           "BPSD": AlleleAggregateObject.BPSD,
                           "MFRAC": AlleleAggregateObject.MFractionAgreed,
                           "MINFRACCALL": AlleleAggregateObject.minFrac,
                           "MINFRACFILTER": minFracAgreedForFilter,
                           "MFA": AlleleAggregateObject.MFA,
                           "MINAAF": minAAF, "MAXAAF": maxAAF,
                           "MQM": AlleleAggregateObject.AveMQ,
                           "MQB": AlleleAggregateObject.AveBQ,
                           "NDP": NDP, "QA": AlleleAggregateObject.SumBQScore,
                           "PFSD": AlleleAggregateObject.PFSD,
                           "NFSD": AlleleAggregateObject.NFSD,
                           "MPF": AlleleAggregateObject.MPF,
                           "TND": AlleleAggregateObject.TND,
                           "MNF": AlleleAggregateObject.MNF,
                           "NDPS": AlleleAggregateObject.NumberDuplexReads,
                           "SHENRANGE": len(flankingBefore),
                           "SHENREF": shenRef,
                           "SHENVAR": shenVar}
        if(TotalCountStr != "default"):
            self.InfoFields["TACS"] = TotalCountStr
        if(TotalFracStr != "default"):
            self.InfoFields["TAFS"] = TotalFracStr
        if(MergedCountStr != "default"):
            self.InfoFields["MACS"] = MergedCountStr
        if(MergedFracStr != "default"):
            self.InfoFields["MAFS"] = MergedFracStr
        if(self.reverseStrandFraction >= 0):
            self.InfoFields["AARSF"] = self.reverseStrandFraction
            # All Alleles Reverse Read Fraction
        if(ampliconFailed >= 0):
            self.InfoFields["NAF"] = ampliconFailed
        self.InfoFields["RSF"] = AlleleAggregateObject.reverseStrandFraction
        self.InfoStr = ";".join(
            [key + "=" + str(self.InfoFields[key])
             for key in sorted(self.InfoFields.keys())])
        self.FormatFields = {"DP": DOC,
                             "DPA": AC,
                             "DPT": DOCTotal,
                             "BQF": FailedBQReads,
                             "AABS": BothStrandAlignment,
                             "MQF": FailedMQReads,
                             "FAF": FailedAFReads,
                             "FQC": FailedQCReads, "FFM": FailedFMReads,
                             "FSR": AlleleAggregateObject.FSR,
                             "AABPSD": AlleleAggregateObject.AABPSD,
                             "AAMBP": AlleleAggregateObject.AAMBP,
                             "MNCS": minNumSS, "MDP": minDuplexPairs,
                             "NDF": FailedNDReads, "MINAF": minAF,
                             "MINFA": AlleleAggregateObject.minFA,
                             "MAXND": AlleleAggregateObject.maxND,
                             "MMQ": AlleleAggregateObject.minMQ,
                             "MBQ": AlleleAggregateObject.minBQ,
                             "MAXNF": AlleleAggregateObject.maxNF,
                             "PVC": MaxPValue,
                             "CONS": self.CONS,
                             "TYPE": "snp",
                             "NUMALL": AlleleAggregateObject.NUMALT}
        ffkeys = sorted(self.FormatFields.keys())
        self.FormatStr = (
            ":".join(ffkeys) +
            "\t" + ":".join(str(
                self.FormatFields[key]) for key in ffkeys))
        self.str = "\t".join(map(str, [self.CHROM,
                                       self.POS, self.ID,
                                       self.CONS, self.ALT,
                                       self.QUAL, self.FILTER,
                                       self.InfoStr, self.FormatStr]))

    def update(self):
        ffkeys = sorted(self.FormatFields.keys())
        self.FormatKey = ":".join(ffkeys)
        self.FormatValue = ":".join([str(self.FormatFields[key])
                                     for key in
                                     ffkeys])
        self.FormatStr = self.FormatKey + "\t" + self.FormatValue
        self.InfoStr = ";".join([key + "=" + str(self.InfoFields[key])
                                 for key in sorted(self.InfoFields.keys())])

    def __str__(self):
        self.update()
        self.str = "\t".join(map(str, [self.CHROM, self.POS,
                                       self.ID, self.REF, self.ALT,
                                       self.QUAL, self.FILTER, self.InfoStr,
                                       self.FormatStr]))
        return self.str
    """

IDInfoDict = {}
IDInfoDict["LEN"] = HeaderInfoLine(
    ID="LEN", Description="Length of allele (deletion or insertion)",
    Type="Integer", Number="A")


IDFilterDict = {}
IDFilterDict["PASS"] = HeaderFilterLine(ID="PASS",
                                        Description="All filters passed")
IDFilterDict["LowComplexity"] = HeaderFilterLine(
    ID="LowComplexity",
    Description=("Variant's flanking regions have Shannon entropy "
                 "below a required threshold."))
IDFilterDict["InsufficientReadPairs"] = HeaderFilterLine(
    ID="InsufficientReadPairs",
    Description=("Variant not supported by sufficient "
                 "concordant pairs of duplex reads"))
IDFormatDict = {}
IDFormatDict["TYPE"] = HeaderFormatLine(
    Type="String", ID="TYPE",
    Description=("The type of allele, either snp, mnp, "
                 "ins, del, tra, or complex."),
    Number="A")
IDFormatDict["NDPS"] = HeaderFormatLine(
    Type="Integer", ID="NDPS", Number="A",
    Description=("Number of duplex reads (where read 1 and read 2 map to the "
                 "same location) supporting variant from both directions."))
IDFormatDict["MDP"] = HeaderFormatLine(
    ID="MDP", Number="A", Type="Float",
    Description=("Mean depth of coverage for region in which indel was "
                 "called. DOC is calculated by averaging 5 bases before and"
                 " 5 bases after the event to give context for the event."))
IDFormatDict["DPA"] = HeaderFormatLine(
    Type="Integer", ID="AC", Number="A",
    Description=("Number of reads passing filters supporting variant."))


def GetIDVCFHeader(fileFormat="default", commandStr="default",
                   reference="default",
                   header="default", sampleName="DefaultSampleName"):
    reference = reference.split("/")[-1]
    # If the reference is the path 
    HeaderLinesStr = ""
    # fileformat line
    HeaderLinesStr += str(HeaderFileFormatLine(
        fileformat=fileFormat)) + "\n"
    # FILTER lines
    for key in sorted(IDFilterDict.keys()):
        HeaderLinesStr += str(IDFilterDict[key]) + "\n"
    # INFO lines
    HeaderLinesStr += "\n".join([str(IDInfoDict[key]) for
                       key in sorted(IDInfoDict.keys())]) + "\n"
    # FORMAT lines
    HeaderLinesStr += "\n".join([str(IDFormatDict[key]) for
                       key in IDFormatDict.keys()]) + "n"
    # commandline line
    if(commandStr != "default"):
        HeaderLinesStr += str(HeaderCommandLine(
            commandStr=commandStr)) + "\n"
    # reference line
    if(reference != "default"):
        HeaderLinesStr += str(HeaderReferenceLine(
            reference=reference)) + "\n"
    # contig lines
    if(header != "default"):
        HeaderLinesStr += GetContigHeaderLines(header) + "\n"
    HeaderLinesStr += "\t".join(["#CHROM", "POS",
                                 "ID", "REF",
                                 "ALT", "QUAL",
                                 "FILTER", "INFO", "FORMAT",
                                 sampleName]) + "\n"
    return HeaderLinesStr