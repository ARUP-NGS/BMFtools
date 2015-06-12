# cython: c_string_type=str, c_string_encoding=ascii

"""
Contains tools related to analyzing tumor/normal pairs.
"""
import pysam
from .Probability import (ConfidenceIntervalAAF, defaultPValue,
                          GetCeiling)
from utilBMF.HTSUtils import makeinfodict, makeformatdict, pl, ThisIsMadness
from subprocess import check_output
import sys
import logging
from math import sqrt as msqrt
asVCF = pysam.asVCF()
cimport pysam.TabProxies
cimport numpy as np
cimport cython
from utilBMF.HTSUtils cimport cystr


@cython.returns(np.longdouble_t)
def GetObservedAAFCeiling(int AC, int DOC,
                          np.longdouble_t pVal=defaultPValue):
    """
    Returns the maximum AAF given observed AC and DOC.
    """
    return ConfidenceIntervalAAF(AC, DOC, pVal=pVal)[1] / DOC


@cython.returns(cython.bint)
def TestTumorVsNormal(int tAC, int tDOC,
                      int nAC, int nDOC,
                      np.longdouble_t pVal=defaultPValue):
    """
    Tests whether the quantity of a variant in a tumor sample analyzed for
    somatic variants is significantly different from the corresponding
    variant in the matched normal.
    PartialPVal takes the 1 / (sqrt(1-p)) so that when it is squared, it
    becomes the pVal of TestTumorVsNormal.
    Handles the numeric
    """
    cdef np.longdouble_t nMaxAAF, tumorMaxNoise
    PartialPVal = 1 - msqrt(1 - pVal)
    nMaxAAF = GetObservedAAFCeiling(nAC, nDOC, pVal=PartialPVal)
    tumorMaxNoise = GetCeiling(tDOC, p=nMaxAAF, pVal=PartialPVal)
    if(tAC / (1. * tDOC) > (tumorMaxNoise)):
        return True
    return False


@cython.returns(cython.bint)
def TumorVNormalVCFProxy(pysam.TabProxies.VCFProxy tumor,
                         pysam.TabProxies.VCFProxy normal=None):
    """
    Tests if a variant in a tumor should be filtered out or not given
    evidence from a paired normal.
    """
    if normal is None:
        raise ThisIsMadness("normal for TumorVNormalVCFProxy must be a "
                            "VCF Proxy Object, not None!")
    tInfoDict = makeinfodict(tumor)
    tFormatDict = makeformatdict(tumor)
    nInfoDict = makeinfodict(normal)
    nFormatDict = makeformatdict(normal)
    return TestTumorVsNormal(tInfoDict["AC"], tFormatDict["DP"],
                             nInfoDict["AC"], nFormatDict["DP"])


def FilterTumorCallsByNormalAAF(tumor, normal="default", outVCF="default",
                                np.longdouble_t pVal=defaultPValue):
    """
    Adds a filter to each record in a tumor VCF with a matched normal.
    If the normal is barcoded, use this pipeline with the same parameters.
    If not, best practice is calling freebayes with
    "--min-alternate-fraction 0" and --pooled-continuous.
    (Will require a different function, actually.)
    """
    cdef pysam.TabProxies.VCFProxy rec, nRec, i
    cdef int nDOC, nAC, nAllelesAtPos
    cdef np.longdouble_t nMaxAAF
    cdef cystr f
    if(normal == "default"):
        raise ThisIsMadness("Noraml VCF required for T/N filtering.")
    if(outVCF == "default"):
        outHandle = sys.stdout
    else:
        outHandle = open(outVCF, "w")
    header = check_output("head -n 4000 %s | grep ^#" % tumor)
    ohw = outHandle.write
    ohw(header)
    tumorIterator = pysam.tabix_iterator(open(tumor, "rb"), asVCF)
    normalHandle = pysam.TabixFile(normal, parser=asVCF)
    nhf = normalHandle.fetch
    for rec in tumorIterator:
        # Load all records with precisely our ref record's position
        try:
            normalRecs = list(nhf(rec.contig + ":" + str(rec.pos - 3) +
                                  "-" + str(rec.pos + 3)))
        except ValueError:
            pl("Looks like contig %s just isn't in the tabix'd" % rec.contig +
               "file. Give up - continuing!", level=logging.DEBUG)
            outVCF.write(str(rec) + "\n")
            continue
        normalRecs = [i for i in normalRecs if i.ref == rec.ref and
                      i.pos == rec.pos]
        nAllelesAtPos = len(normalRecs)
        normalRecs = [i for i in normalRecs if i.alt == rec.alt]
        # Get just the record (or no record) that has that alt.
        if(nAllelesAtPos == 0):
            pl("No variants called at normal position. %s" % str(rec),
               level=logging.DEBUG)
            outVCF.write(str(rec) + "\n")
            continue
        if(len(normalRecs) == 0):
            pl("Looks like this variant: "
               "%s wasn't called at all in the normal." % (str(rec)),
               level=logging.DEBUG)
            outVCF.write(str(rec) + "\n")
            continue
        nRec = normalRecs[0]
        VariantPasses = TestTumorVsNormal(rec, normal=nRec)
        if(not VariantPasses):
            rFilterList = rec.filter.split(";")
            if("PASS" in rFilterList):
                rec.filter = ";".join(sorted([f for f in rFilterList if
                                              f != "PASS"] +
                                             ["TumorNormalIndistinguishable"]))
        outVCF.write(str(rec) + "\n")
    return outVCF
