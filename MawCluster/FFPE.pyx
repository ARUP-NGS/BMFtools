# cython: c_string_type=str, c_string_encoding=ascii
# cython: profile=True, cdivision=True, cdivision_warnings=True

from cytoolz import map as cmap, memoize
from math import log10 as mlog10
from .BCVCF import IterativeVCFFile
from .Probability import ConfidenceIntervalAAF, GetCeiling
from .SNVUtils import HeaderFilterDict, HeaderFunctionCallLine
from operator import methodcaller as mc
from subprocess import check_output
from utilBMF.ErrorHandling import IllegalArgumentError, ThisIsMadness
import cStringIO
import cython
import logging
import math
import numpy as np
import operator
import pysam
import shlex
import subprocess
import uuid
from utilBMF.HTSUtils import (printlog as pl,
                              NameSortAndFixMate, makeinfodict,
                              MakeVCFProxyDeaminationFilter)
from BMFMain.main import __version__ as BMFVersion
cimport pysam.TabProxies
cimport numpy as np
cimport cython


"""
Contains utilities relating to FFPE
"""


@cython.locals(maxFreq=np.longdouble_t)
@cython.returns(np.longdouble_t)
def GetDeaminationFrequencies(inVCF, maxFreq=0.15, FILTER="",
                              minGCcount=50):

    """
    Returns a list of raw base frequencies for G->A and C->T.
    Only accepts SNVCrawler's VCF output as it requires those INFO fields.
    FILTER must be a comma-separated list of FILTER strings as defined
    in the VCF header.
    This is slower than the Py version, so I will be working with it to
    get some speed soon.
    """
    cdef cython.int TotalCG
    cdef cython.int TotalCG_TA
    cdef cython.int DP
    cdef cython.int GC
    cdef np.longdouble_t freq
    validFilters = list(cmap(mc("lower"), HeaderFilterDict.keys()))
    FILTER = FILTER.lower()
    filters = FILTER.split(",")
    pl("Filters: %s" % filters)
    pl("maxFreq: %s" % maxFreq)
    if(FILTER != ""):
        for filter in filters:
            if filter not in validFilters:
                raise IllegalArgumentError(
                    "Filter must be a valid VCF Filter."
                    "Valid filters: %s" % validFilters)
    IVCFObj = IterativeVCFFile(inVCF)
    TotalCG = 0
    TotalCG_TA = 0
    for line in IVCFObj:
        if(FILTER != ""):
            for filter in filters:
                if filter in line.FILTER.lower():
                    pl("Failing line for filter %s" % filter,
                       level=logging.DEBUG)
                    continue
        lid = line.InfoDict
        cons = lid["CONS"]
        MACSStr = lid["MACS"]
        macsDict = dict(list(cmap(mc(
            "split", ">"), MACSStr.split(","))))
        if(cons == "C" and line.REF == "C"):
            GC = int(macsDict["C"])
            TA = int(macsDict["T"])
            if(GC < minGCcount):
                continue
            if(GC == 0):
                continue
            if(TA * 1.0 / GC <= maxFreq):
                TotalCG_TA += TA
                TotalCG += GC
        elif(cons == "G" or line.REF == "G"):
            GC = int(macsDict["G"])
            if(GC < minGCcount):
                continue
            TA = int(macsDict["A"])
            if(GC == 0):
                continue
            if(TA * 1.0 / GC <= maxFreq):
                TotalCG_TA += TA
                TotalCG += GC
    freq = 1. * TotalCG_TA / TotalCG
    pl("TotalCG_TA: %s. TotalCG: %s." % (TotalCG_TA, TotalCG))
    pl("Estimated frequency: %s" % freq)
    pl("For perspective, a 0.001 pValue ceiling at 100 DOC would be "
       "%s" % (GetCeiling(100, p=freq, pVal=0.001) / 100.))
    pl("Whereas, a 0.001 pValue ceiling at 1000 DOC would be %s" % (GetCeiling(
        1000, p=freq, pVal=0.001) / 1000.))
    return freq


@cython.returns(np.longdouble_t)
def PyGetDeamFreq(inVCF, maxFreq=0.15, FILTER="",
                  minGCcount=50):

    """
    Trying to see if the C compilation is throwing anything off.
    Returns a list of raw base frequencies for G->A and C->T.
    Only accepts SNVCrawler's VCF output as it requires those INFO fields.
    FILTER must be a comma-separated list of FILTER strings as defined
    in the VCF header.
    """
    validFilters = list(cmap(mc("lower"), HeaderFilterDict.keys()))
    filters = FILTER.lower().split(",")
    pl("Filters: %s" % filters)
    pl("maxFreq: %s" % maxFreq)
    if(FILTER != ""):
        for filter in filters:
            if filter not in validFilters:
                raise IllegalArgumentError(
                    "Filter must be a valid VCF Filter."
                    " %s" % validFilters)
    IVCFObj = IterativeVCFFile(inVCF)
    TotalCG = 0
    TotalCG_TA = 0
    for line in IVCFObj:
        if(FILTER != ""):
            for filter in filters:
                if filter in line.FILTER.lower():
                    pl("Failing line for filter %s" % filter,
                       level=logging.DEBUG)
                    continue
        lid = line.InfoDict
        cons = lid["CONS"]
        MACSStr = lid["MACS"]
        macsDict = dict(list(cmap(mc(
            "split", ">"), MACSStr.split(","))))
        if((cons == "C" or line.REF == "C") and line.ALT == "C"):
            GC = int(macsDict["C"])
            TA = int(macsDict["T"])
            if(GC < minGCcount):
                continue
            if(GC == 0):
                continue
            if(TA * 1.0 / GC <= maxFreq):
                TotalCG_TA += TA
                TotalCG += GC
        elif((cons == "G" or line.REF == "G") and line.ALT == "A"):
            GC = int(macsDict["G"])
            if(GC < minGCcount):
                continue
            TA = int(macsDict["A"])
            if(GC == 0):
                continue
            if(TA * 1.0 / GC <= maxFreq):
                TotalCG_TA += TA
                TotalCG += GC
    freq = 1. * TotalCG_TA / TotalCG
    pl("TotalCG_TA: %s. TotalCG: %s." % (TotalCG_TA, TotalCG))
    pl("Estimated frequency: %s" % freq)
    pl("For perspective, a 0.001 pValue ceiling at 100 DOC would be "
       "%s" % (GetCeiling(100, p=freq, pVal=0.001) / 100.))
    pl("Whereas, a 0.001 pValue ceiling at 1000 DOC would be %s" % (GetCeiling(
        1000, p=freq, pVal=0.001) / 1000.))
    return freq


@cython.locals(pVal=np.longdouble_t, DOC=cython.int,
               maxFreqNoise=np.longdouble_t, ctfreq=np.longdouble_t,
               AAF=np.longdouble_t, recordsPerWrite=cython.int)
def FilterByDeaminationFreq(inVCF, pVal=0.001, ctfreq=0.018,
                            recordsPerWrite=5000, outVCF="default"):
    """
    If observed AAF is greater than the upper limit of the confidence window
    with a given P-Value, the variant is permitted to stay.
    Otherwise, DeaminationNoise replaces PASS or is appended to other filters.
    """
    pl("C-T/G-A frequency set to %s" % ctfreq)
    IVCFObj = IterativeVCFFile(inVCF)
    if(outVCF == "default"):
        outVCF = ".".join(inVCF.split(".")[0:-1] + ["ctfilt", "vcf"])
    pl("FilterByDeaminationFreq called. inVCF: %s. outVCF: %s." % (inVCF,
                                                                   outVCF))
    outHandle = open(outVCF, "w")
    mfdnpStr = str(int(-10 * mlog10(pVal)))
    functionCall = ("FilterByDeaminationFreq(%s, pVal=%s, " % (inVCF, pVal) +
                    "ctfreq=%s, recordsPerWrite=" % ctfreq +
                    "%s). BMFTools version: %s" % (recordsPerWrite,
                                                   BMFVersion))
    IVCFObj.header.insert(-1, str(HeaderFunctionCallLine(functionCall)))
    outHandle.write("\n".join(IVCFObj.header) + "\n")
    recordsArray = []
    for line in IVCFObj:
        if(len(recordsArray) >= recordsPerWrite):
            outHandle.write("\n".join(map(str, recordsArray)) +
                            "\n")
            recordsArray = []
        if(line.REF != "C" or line.REF != "G"):
            recordsArray.append(line)
            continue
        if(line.REF == "C" and line.ALT != "T"):
            recordsArray.append(line)
            continue
        if(line.REF == "G" and line.ALT != "A"):
            recordsArray.append(line)
            continue
        DOC = int(line.GenotypeDict["DP"])
        AAF = float(line.InfoDict["AC"]) / DOC
        maxFreqNoise = GetCeiling(DOC, p=ctfreq, pVal=pVal) / (DOC * 1.)
        if AAF < maxFreqNoise:
            if(line.FILTER == "PASS"):
                line.FILTER = "DeaminationNoise"
            else:
                line.FILTER += ",DeaminationNoise"
        line.InfoDict["MFDN"] = maxFreqNoise
        line.InfoDict["MFDNP"] = mfdnpStr
        recordsArray.append(line)
    outHandle.write("\n".join(map(str, recordsArray)) + "\n")
    outHandle.flush()
    outHandle.close()
    return outVCF


@cython.locals(maxFreq=np.longdouble_t, pVal=np.longdouble_t)
def TrainAndFilter(inVCF, maxFreq=0.1, FILTER="",
                   pVal=0.001):
    """
    Calls both GetDeaminationFrequencies and FilterByDeaminationFreq.
    """
    cdef np.longdouble_t DeamFreq
    DeamFreq = GetDeaminationFrequencies(inVCF, maxFreq=maxFreq,
                                         FILTER=FILTER)
    pl("Estimated deamination frequency: %s" % DeamFreq)
    OutVCF = FilterByDeaminationFreq(inVCF, pVal=pVal,
                                     ctfreq=DeamFreq)
    pl("Output VCF: %s" % OutVCF)
    return OutVCF


@cython.locals(primerLen=cython.int, fixmate=cython.bint)
def PrefilterAmpliconSequencing(inBAM, primerLen=20, outBAM="default",
                                fixmate=True):
    """
    This program outputs a BAM file which eliminates potential mispriming
    events from the input BAM file.
    """
    if(outBAM == "default"):
        outBAM = ".".join(inBAM.split(".")[0:-1] + ["amplicon",
                                                    "filt", "bam"])
    pl("Primer length set to %s for prefiltering." % primerLen)
    pl("OutBAM: %s" % outBAM)
    pl("fixmate: %s" % fixmate)
    tempFile = str(uuid.uuid4().get_hex().upper()[0:8]) + ".bam"
    inHandle = pysam.AlignmentFile(inBAM, "rb")
    outHandle = pysam.AlignmentFile(tempFile, "wb", template=inHandle)
    for rec in inHandle:
        tempQual = rec.qual[primerLen:]
        rec.seq = rec.seq[primerLen:]
        rec.qual = tempQual
        if(rec.is_reverse):
            rec.pos -= primerLen
        else:
            rec.pos += primerLen
        outHandle.write(rec)
    inHandle.close()
    outHandle.close()
    if(fixmate):
        newTemp = NameSortAndFixMate(tempFile, sortAndIndex=True)
        subprocess.check_call(["mv", newTemp, outBAM])
    return outBAM


@memoize
@cython.locals(l=cython.str)
@cython.returns(cython.float)
def getFreq(pysam.TabProxies.VCFProxy rec, l="d"):
    """
    Returns allele frequency for a tabix VCF Proxy made from SNVCrawler.
    """
    return float(dict([i.split(">") for i in
                       makeinfodict(rec)["MAFS"].split(",")])[l])


@cython.returns(cython.float)
def GetTabixDeamFreq(cython.str inVCF):
    """
    Gets deamination frequency for a tabixed VCF file, under the assumption
    that the majority of C-T/G-A calls at low frequencies which are not
    ablated by family demultiplexing are due to formalin fixation.
    """
    cdef cython.int atCounts
    cdef cython.int gcCounts
    cdef pysam.TabProxies.VCFProxy rec
    cdef cython.float freq
    atCounts = 0
    gcCounts = 0
    import sys
    import pysam
    a = pysam.tabix_iterator(open(inVCF, "rb"), pysam.asVCF())
    for rec in a:
        mid = makeinfodict(rec)
        if(mid["CONS"] == "C" and
           getFreq(rec, "T") / getFreq(rec, "C") < 0.15 and
           rec.alt == "T"):
            atCounts += int(dict([i.split(">") for
                                  i in mid["MACS"].split(",")])["T"])
            gcCounts += int(dict([i.split(">") for
                                  i in mid["MACS"].split(",")])["C"])
        if(mid["CONS"] == "G" and
           getFreq(rec, "A") / getFreq(rec, "G") < 0.15 and
           rec.alt == "A"):
            atCounts += int(dict([i.split(">") for
                                  i in mid["MACS"].split(",")])["A"])
            gcCounts += int(dict([i.split(">") for
                                  i in mid["MACS"].split(",")])["G"])
        if(rec.ref == "C" and getFreq(rec, "T") < 0.25 and
           getFreq(rec, "C") >= 0.3 and rec.alt == "T"):
            atCounts += int(dict([i.split(">") for
                                  i in mid["MACS"].split(",")])["T"])
            gcCounts += int(dict([i.split(">") for
                                  i in mid["MACS"].split(",")])["C"])
        if(rec.ref == "G" and getFreq(rec, "A") < 0.25 and
           getFreq(rec, "G") >= 0.3 and rec.alt == "A"):
            atCounts += int(dict([i.split(">") for
                                  i in mid["MACS"].split(",")])["A"])
            gcCounts += int(dict([i.split(">") for
                                  i in mid["MACS"].split(",")])["G"])
    freq = (1. * atCounts) / gcCounts
    print("Final atCounts: %s" % atCounts)
    print("Final gcCounts: %s" % gcCounts)
    print("Est deam freq: %s" % (freq))
    return freq


@cython.locals(pVal=np.longdouble_t, DOC=cython.int,
               maxFreqNoise=np.longdouble_t, ctfreq=np.longdouble_t,
               AAF=np.longdouble_t, recordsPerWrite=cython.int)
def TabixDeamFilter(inVCF, pVal=0.001, ctfreq=0.006,
                    recordsPerWrite=5000, outVCF="default"):
    """
    If observed AAF is greater than the upper limit of the confidence window
    with a given P-Value, the variant is permitted to stay.
    Otherwise, DeaminationNoise replaces PASS or is appended to other filters.
    """
    pl("C-T/G-A frequency set to %s" % ctfreq)
    inHandle = pysam.tabix_iterator(open(inVCF, "rb"), pysam.asVCF())
    if(outVCF == "default"):
        outVCF = ".".join(inVCF.split(".")[0:-2] + ["ctfilt", "vcf"])
    headerStringIO = cStringIO.StringIO()
    headerStringIO.write(check_output("zcat %s | head -n 2000" % inVCF,
                                      shell=True))
    headerStringIO.reset()
    headerLines = IterativeVCFFile(headerStringIO).header
    del headerStringIO
    pl("TabixDeamFilter called. inVCF: %s. outVCF: %s." % (inVCF,
                                                           outVCF))
    if(not isinstance(outVCF, file)):
        outHandle = open(outVCF, "w")
    else:
        outHandle = outVCF
    ohw = outHandle.write
    mfdnpStr = str(int(-10 * mlog10(pVal)))
    functionCall = ("FilterByDeaminationFreq(%s, pVal=%s, " % (inVCF, pVal) +
                    "ctfreq=%s, recordsPerWrite=" % ctfreq +
                    "%s). BMFTools version: %s" % (recordsPerWrite,
                                                   BMFVersion))
    headerLines.insert(-1, str(HeaderFunctionCallLine(functionCall)))
    ohw("\n".join(headerLines) + "\n")
    FilterFn = MakeVCFProxyDeaminationFilter(ctfreq, conf=pVal,
                                             key="MFDNP",
                                             value=mfdnpStr)
    recordsArray = []
    for rec in inHandle:
        recordsArray.append(FilterFn(rec))
        if(len(recordsArray) >= recordsPerWrite):
            outHandle.write("\n".join(list(cmap(str, recordsArray))) + "\n")
            recordsArray = []
    outHandle.write("\n".join(list(cmap(str, recordsArray))) + "\n")
    outHandle.flush()
    outHandle.close()
    return outVCF
