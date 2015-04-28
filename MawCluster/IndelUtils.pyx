# cython: c_string_type=str, c_string_encoding=ascii
# cython: profile=True, cdivision=True, cdivision_warnings=True
import pysam
import cython
from cytoolz import map as cmap
from utilBMF.HTSUtils import (FractionAligned, FractionSoftClipped,
                              printlog as pl)
from utilBMF.ErrorHandling import ThisIsMadness
cimport pysam.calignmentfile

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
    if(indelOutputBAM == "default"):
        indelOutputBAM = ".".join(inBAM.split(".")[0:-1] + ["indRel", "bam"])
    if(otherOutputBAM == "default"):
        otherOutputBAM = ".".join(inBAM.split(".")[0:-1] + ["indIrl", "bam"])
    inHandle = pysam.AlignmentFile(inBAM, "rb")
    indelHandle = pysam.AlignmentFile(indelOutputBAM, "wb", template=inHandle)
    otherHandle = pysam.AlignmentFile(otherOutputBAM, "wb", template=inHandle)
    ohw = otherHandle.write
    ihw = indelHandle.write
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
        else:
            ohw(read1)
            ohw(read2)
    inHandle.close()
    otherHandle.close()
    indelHandle.close()
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
    try:
        if(read.opt("SV") != "NF"):
            return True
    except KeyError:
        # No SV tag was found, not much we can do.
        pass
    return False


def GetSCFractionArray(inBAM):
    cdef pysam.calignmentfile.AlignedSegment i
    inHandle = pysam.AlignmentFile(inBAM, "rb")
    return [FractionSoftClipped(i) for i in inHandle]


def GetFreebayesCallStr(inBAM, ref="default", bed="default",
                        outVCF="default", ploidy=-1, K=True,
                        minMQ=1, minBQ=3, haplotypeLength=80,
                        rdf=1.0, bestNAlleles=10):
    """
    Used to call freebayes for indel calling.
    """
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
                "--min-alternate-fraction 0 --pooled-continuous")
    else:
        pl("Ploidy set: %s" % ploidy)
        cStr = ("freebayes -v %s -t %s -f %s " % (outVCF, bed, ref) +
                "-q %s -m %s --haplotype-length" % (minBQ, minMQ) +
                " %s -D %s " % (haplotypeLength, rdf) +
                "--min-alternate-fraction 0 --pooled-continuous" +
                " -p %s --use-best-n-alleles %s" % (ploidy, bestNAlleles))
    return cStr


@cython.returns(cython.str)
def GetFBOutVCFFromStr(cython.str cStr):
    """
    Gets out vcf from freebayes call. Used for parallelization.
    """
    return cStr.split(" ")[2]
