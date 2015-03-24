# cython: c_string_type=str, c_string_encoding=ascii
# cython: profile=True, cdivision=True, cdivision_warnings=True

import logging
import operator
import math
from math import log10 as mlog10
from operator import methodcaller as mc

import cython
cimport cython
import numpy as np
cimport numpy as np

from MawCluster.BCVCF import IterativeVCFFile
from utilBMF.HTSUtils import printlog as pl, ThisIsMadness
from utilBMF.ErrorHandling import IllegalArgumentError
from MawCluster.SNVUtils import HeaderFilterDict, HeaderFunctionCallLine
from MawCluster.Probability import ConfidenceIntervalAAF, GetCeiling
from MawCluster.BCVCF import IterativeVCFFile

ctypedef np.longdouble_t dtype128_t


"""
Contains utilities relating to FFPE
"""

BMFVersion = "0.0.7.2"


@cython.locals(pVal=dtype128_t, DOC=cython.long,
               maxFreqNoise=dtype128_t, ctfreq=dtype128_t, AAF=dtype128_t,
               recordsPerWrite=cython.long)
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
    IVCFObj.header.insert(-1, HeaderFunctionCallLine(functionCall).ToString())
    outHandle.write("\n".join(IVCFObj.header) + "\n")
    recordsArray = []
    for line in IVCFObj:
        if(len(recordsArray) >= recordsPerWrite):
            outHandle.write("\n".join(map(mc("ToString"), recordsArray)) +
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
    outHandle.write("\n".join(map(mc("ToString"), recordsArray)) + "\n")
    outHandle.flush()
    outHandle.close()
    return outVCF


@cython.locals(maxFreq=dtype128_t)
def GetDeaminationFrequencies(inVCF, maxFreq=0.2,
                              FILTER=""):

    """
    Returns a list of raw base frequencies for G->A and C->T.
    Only accepts SNVCrawler's VCF output as it requires those INFO fields.
    FILTER must be a comma-separated list of FILTER strings as defined
    in the VCF header.
    """
    cdef cython.long TotalCG
    cdef cython.long TotalCG_TA
    cdef cython.long DP
    cdef cython.long GC
    validFilters = map(mc("lower"), HeaderFilterDict.keys())
    FILTER = FILTER.lower()
    filters = FILTER.split(",")
    pl("Filters: %s" % filters)
    pl("maxFreq: %s" % maxFreq)
    if(FILTER != ""):
        for filter in filters:
            if filter not in validFilters:
                raise IllegalArgumentError("Filter must be a valid VCF Filter. "
                                           "%s" % validFilters)
    IVCFObj = IterativeVCFFile(inVCF)
    TotalCG = 0
    TotalCG_TA = 0
    for line in IVCFObj:
        if(FILTER != ""):
            for filter in filters:
                if filter in line.FILTER.lower():
                    pl("Failing line for filter %s" % filter, level=logging.DEBUG)
                    continue
        if(line.REF == "C"):
            if(line.ALT != "T"):
                continue
            GC = int(line.InfoDict["MACS"].split(",")[1].split(">")[1])
            TA = int(line.InfoDict["MACS"].split(",")[2].split(">")[1])
            if(GC == 0):
                continue
            if(TA * 1.0 / GC <= maxFreq):
                TotalCG_TA += TA
                TotalCG += GC
        elif(line.REF == "G"):
            if(line.ALT != "A"):
                continue
            GC = int(line.InfoDict["MACS"].split(",")[3].split(">")[1])
            TA = int(line.InfoDict["MACS"].split(",")[0].split(">")[1])
            if(GC == 0):
                continue
            if(TA * 1.0 / GC <= maxFreq):
                TotalCG_TA += TA
                TotalCG += GC
    return 1. * TotalCG_TA / TotalCG



@cython.locals(maxFreq=dtype128_t, pVal=dtype128_t)
def TrainAndFilter(inVCF, maxFreq=0.1, FILTER="",
                   pVal=0.001):
    """
    Calls both GetDeaminationFrequencies and FilterByDeaminationFreq.
    """
    cdef np.ndarray[dtype128_t, ndim = 1] DeamFreqs
    cdef dtype128_t DeamFreq
    DeamFreqs = GetDeaminationFrequencies(inVCF, maxFreq=maxFreq,
                                         FILTER=FILTER)
    DeamFreq = np.mean(DeamFreqs, dtype=np.longdouble)
    pl("Estimated deamination frequency: %s" % DeamFreq)
    OutVCF = FilterByDeaminationFreq(inVCF, pVal=pVal,
                                     ctfreq=DeamFreq)
    pl("Output VCF: %s" %OutVCF)
    return OutVCF
