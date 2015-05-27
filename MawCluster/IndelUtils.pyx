# cython: c_string_type=str, c_string_encoding=ascii
# cython: profile=True, cdivision=True, cdivision_warnings=True
import operator
import pysam
import cython
from cytoolz import map as cmap
from utilBMF.HTSUtils import (FractionAligned, FractionSoftClipped,
                              printlog as pl, CoorSortAndIndexBam,
                              GetKmersToCheck, FastqStrFromKmerList,
                              Bowtie2FqToStr, GetMQPassReads,
                              GetInsertionFromAlignedSegment,
                              GetDeletionFromAlignedSegment,
                              FacePalm, ParseBed,
                              shen, ssStringFromRead, ccopy,
                              Insertion, Deletion, IndelQuiver,
                              AbstractIndelContainer, IDVCFLine)
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
ctypedef utilBMF.HTSUtils.Insertion Insertion_t
ctypedef utilBMF.HTSUtils.Deletion Deletion_t

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
    cdef cython.int idrel
    cdef cython.int idirl
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


@cython.returns(cython.bint)
def IsIndelRelevant(
        pysam.calignmentfile.AlignedSegment read, cython.int minFam=2,
        cython.float minSF=0.6, cython.bint keepUnmapped=False):
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
def GetUniquelyMappableKmers(cython.str ref, cython.int k=30,
                             list bedline=[], cython.int minMQ=1,
                             cython.int padding=-1,
                             cython.int mismatches=-1):
    """
    Uses a set of HTSUtils methods to find kmers from a region
    which are uniquely mappable. This makes it possible to do alignment-free
    variant-calling. (Well, except for the bwasw step).
    If no outfile is specified, defaults to stdout
    """
    cdef list kmerList, PassingReadNames
    cdef cython.str fqStr, bowtieStr
    pl("Getting potential kmers for k=%s" % k)
    kmerList = GetKmersToCheck(ref, k=k, bedline=bedline, padding=padding)
    pl("Making dummy fastq records for each kmer")
    fqStr = FastqStrFromKmerList(kmerList)
    pl("Aligning these kmers to the genome to test for unique mappability"
       " with a given number minMQ %s." % minMQ)
    bowtieStr = Bowtie2FqToStr(fqStr, ref=ref, seed=k, mismatches=mismatches)
    PassingReadNames = GetMQPassReads(bowtieStr, minMQ=minMQ)
    return PassingReadNames


@cython.returns(IndelQuiver_t)
def FillIndelQuiverRegion(inBAM, cython.int minPairs=2,
                          cython.float minShen=0.2,
                          cython.int window=16, cython.str ref=None,
                          list bedRegion=[], cython.int minMQ=1,
                          cython.int minFM=1, cython.int minNumSS=1):
    """
    To fill an IndelQuiver, the bam must have the SV tags set, and positive
    calls require a DSI or DSD tag.

    You do need to be somewhat careful about adding Insertion or Deletion
    objects to an IndelQuiver object, as those objects could have been made
    with reads not satisfying minMQ and minFM.
    """
    cdef IndelQuiver_t Quiver
    cdef pysam.calignmentfile.AlignedSegment rec
    cdef pysam.calignmentfile.AlignmentFile inHandle
    cdef cython.int end
    Quiver = IndelQuiver(ref=ref, window=window, minFM=minFM, minMQ=minMQ,
                         bam=inBAM, minShen=minShen, minPairs=minPairs,
                         minNumSS=minNumSS)
    refHandle = pysam.FastaFile(ref)
    inHandle = pysam.AlignmentFile(inBAM, "rb")
    inFetch = inHandle.fetch(bedRegion[0], bedRegion[1])
    end = bedRegion[2]
    for rec in inFetch:
        if(rec.mapping_quality < minMQ):
            continue
            """
            Since we're looking for DS[ID], if one read has a mapq
            sufficient for inclusion, we get the benefit of seeing that indel
            even if one of the reads failed the MQ threshold.
            """
        if(rec.opt("FM") < minFM):
            continue
            """
            Only working with read families with FM >= minFM
            cleans up our results.
            """
        try:
            svtags = rec.opt("SV")
        except KeyError:
            FacePalm("SV tags must be marked to call tagged indels!")
        if("DSI" in svtags):
            Quiver.addIndel(GetInsertionFromAlignedSegment(rec,
                                                           handle=refHandle))
        if("DSD" in svtags):
            Quiver.addIndel(GetDeletionFromAlignedSegment(rec,
                                                          handle=refHandle))
        if(rec.reference_pos >= end):
            pl("Finished bed region %s." % bedRegion)
            break
    return Quiver


@cython.returns(IndelQuiver_t)
def FillEntireQuiver(inBAM, cython.int minPairs=2, cython.float minShen=0.2,
                     cython.int window=16, cython.str ref=None,
                     cython.str bed=None, cython.int minFM=1,
                     cython.int minMQ=1, cython.int minNumSS=1):
    """
    Simply fills a quiver for each region, merges the quivers, and returns
    the full quiver.
    """
    cdef IndelQuiver_t Quiver, RegionalQuiver
    cdef list bedlines, line
    Quiver = IndelQuiver(ref=ref, window=window, minFM=minFM, minMQ=minMQ,
                         bam=inBAM, minShen=minShen, minPairs=minPairs,
                         minNumSS=minNumSS)
    bedlines = ParseBed(bed)
    for line in bedlines:
        RegionalQuiver = FillIndelQuiverRegion(inBAM, window=window,
                                               minFM=minFM, minMQ=minMQ,
                                               minPairs=minPairs,
                                               minShen=minShen, ref=ref,
                                               bedRegion=line,
                                               minNumSS=minNumSS)
        Quiver.mergeQuiver(RegionalQuiver)
    return Quiver


# Info fields
IDInfoDict = {}
IDInfoDict["MINSHEN"] = HeaderInfoLine(
    ID="MINSHEN", Type="Float", Number="1",
    Description=("Minimum Shannon entropy for both"
                 " preceding and succeeding regions"))
IDInfoDict["SHENWINDOW"] = HeaderInfoLine(
    ID="SHENWINDOW", Type="Integer", Number="1",
    Description=("Number of preceding or succeeding bases to include in "
                 "Shannon entropy calculations."))


# Filter fields
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
IDFilterDict["InsufficientStopStarts"] = HeaderFilterLine(
    ID="InsufficientStopStarts",
    Description="Number of start/stop coordinates is less than minNumSS")

# Format fields
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
IDFormatDict["DPA"] = HeaderFormatLine(
    Type="Integer", ID="AC", Number="A",
    Description=("Number of reads passing filters supporting variant."))
IDFormatDict["SHEN"] = HeaderFormatLine(
    ID="SHEN", Type="Float", Number="A",
    Description=("Shannon entropy of preceding or succeeding "
                 "sequence, whichever is smaller."))
IDFormatDict["LEN"] = HeaderFormatLine(
    ID="LEN", Description="Length of allele (deletion or insertion)",
    Type="Integer", Number="A")
IDFormatDict["MDP"] = HeaderFormatLine(
    ID="MDP", Number="A", Type="Float",
    Description=("Mean depth of coverage for region in which indel was "
                 "called. DOC is calculated by averaging 5 bases before and"
                 " 5 bases after the event to give context for the event."))


def GetIDVCFHeader(fileFormat="default", commandStr="default",
                   reference="default",
                   header="default", sampleName="DefaultSampleName"):
    reference = reference.split("/")[-1]
    # If the reference is a path, trim to just the name of the file.
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