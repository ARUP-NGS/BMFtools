# cython: c_string_type=str, c_string_encoding=ascii
# cython: profile=True, cdivision=True, cdivision_warnings=True
import pysam
import cython
from cytoolz import map as cmap
from utilBMF.HTSUtils import (FractionAligned, FractionSoftClipped,
                              printlog as pl, CoorSortAndIndexBam,
                              GetKmersToCheck, FastqStrFromKmerList,
                              BowtieFqToStr, GetMQPassReads)
from utilBMF.ErrorHandling import ThisIsMadness
import os.path
import uuid
from cytoolz.itertoolz import frequencies as cyfreq
import sys
cimport pysam.calignmentfile
cimport cython

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
        cython.float minSF=0.2, cython.bint keepUnmapped=False):
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