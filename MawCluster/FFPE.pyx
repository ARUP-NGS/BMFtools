# cython: c_string_type=str, c_string_encoding=ascii
# cython: profile=True, cdivision=True, cdivision_warnings=True

import logging
import operator
import math
from math import log10 as mlog10

import cython
cimport cython
import numpy as np
cimport numpy as np
from sklearn.covariance import EllipticEnvelope
from sklearn import svm
from scipy.misc import comb

from MawCluster.BCVCF import IterativeVCFFile
from utilBMF.HTSUtils import printlog as pl, ThisIsMadness
from utilBMF.ErrorHandling import IllegalArgumentError
from MawCluster.SNVUtils import HeaderFilterDict, HeaderFunctionCallLine
from MawCluster.Probability import ConfidenceInterval_, defaultPValue
from MawCluster.BCVCF import IterativeVCFFile
from BMFMain.main import __version__ as BMFVersion

ctypedef np.longdouble_t dtype128_t


"""
Contains utilities relating to FFPE 
"""


@cython.locals(p=dtype128_t, pVal=dtype128_t,
               n=cython.long)
@cython.returns(dtype128_t)
def GetMax_CT_AG(doc, p=0.0, pVal=defaultPValue):
    """
    Given a depth of coverage, an estimated fixation cytosine deamination
    frequency, and a p value, computes the minimum C-T/G-A frequency
    to not be considered an FFPE artefact.
    p: frequency of deamination
    doc: depth of coverage at position
    """
    return ConfidenceInterval_(doc, p=p, pVal=pVal)[1]


@cython.locals(observedFreq=dtype128_t, pVal=dtype128_t,
               doc=cython.long)
@cython.returns(np.ndarray)
def GetAAFBounds_(doc, observedFreq=0.0, pVal=0.05):
    """
    Given a depth of coverage, an observed frequency, and a p value,
    returns a list containing the minimum and maximum AA frequencies,
    respectively, for the actual AAF, within the given p value.
    doc: depth of coverage at position
    default pVal: 0.05
    """
    return np.array(ConfidenceInterval_(doc, p=observedFreq, pVal=pVal)[:2],
                    dtype=np.longdouble)


@cython.locals(numSD=cython.long, pVal=dtype128_t, DOC=cython.long,
               maxFreqNoise=dtype128_t, ctfreq=dtype128_t, AAF=dtype128_t,
               mfdnp=cython.long, recordsPerWrite=cython.long)
def FilterByDeaminationFreq(inVCF, pVal=0.001, ctfreq=0.018,
                            recordsPerWrite=5000):
    """
    If observed AAF is greater than the upper limit of the confidence window
    with a given P-Value, the variant is permitted to stay.
    Otherwise, DeaminationNoise replaces PASS or is appended to other filters.
    """
    pl("C-T/G-A frequency set to %s" % ctfreq)
    IVCFObj = IterativeVCFFile(inVCF)
    outVCF = ".".join(inVCF.split(".")[0:-1] + ["ctfilt", "vcf"])
    outHandle = open(outVCF, "w")
    mfdnp = int(-10 * mlog10(pVal))
    functionCall = ("FilterByDeaminationFreq(%s, pVal=%s, " % (inVCF, pVal) +
                    "ctfreq=%s, recordsPerWrite=" % ctfreq +
                    "%s). BMFTools version: %s" (recordsPerWrite, BMFVersion))
    IVCFObj.header = IVCFObj.header.insert(len(IVCFObj.header) - 1,
        HeaderFunctionCallLine(functionCall))
    outHandle.write("\n".join(IVCFObj.header) + "\n")
    recordsArray = []
    for line in IVCFObj:
        if(len(recordsArray) >= recordsPerWrite):
            outHandle.write("\n".join([line.ToString() for line in recordsArray]) + "\n")
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
        AAF = float(line.InfoDict["AF"])
        DOC = int(line.GenotypeDict["DP"])
        maxFreqNoise = ConfidenceInterval_(DOC, p=ctfreq, pVal=pVal)[1]
        if AAF < maxFreqNoise:
            if(line.FILTER == "PASS"):
                line.FILTER = "DeaminationNoise"
            else:
                line.FILTER += ",DeaminationNoise"
        line.InfoDict["MFDN"] = maxFreqNoise
        line.InfoDict["MFDNP"] = mfdnp
        recordsArray.append(line)
    outHandle.write("\n".join([line.ToString() for line in recordsArray]) + "\n")
    


@cython.locals(FSD=dtype128_t, maxFreq=dtype128_t,
               concatenate=cython.bint)
def GetDeaminationFrequencies(inVCF, maxFreq=0.05, FILTER="LowQual",
                              concatenate=True):

    """
    Returns a list of raw base frequencies for G->A and C->T.
    Only accepts SNVCrawler's VCF output as it requires those INFO fields.
    FILTER must be a comma-separated list of FILTER strings as defined
    in the VCF header.
    """
    validFilters = [i.lower() for i in HeaderFilterDict.keys()]
    FILTER = FILTER.lower()
    filters = FILTER.split(",")
    for filter in filters:
        if filter not in validFilters:
            raise IllegalArgumentError("Filter must be a valid VCF Filter. "
                                       "%s" % validFilters)
    cdef np.ndarray[dtype128_t, ndim = 1] GAFreqNP
    cdef np.ndarray[dtype128_t, ndim = 1] CTFreqNP
    IVCFObj = IterativeVCFFile(inVCF)
    GAFreqArray = []
    CTFreqArray = []
    for line in IVCFObj:
        for filter in filters:
            if filter in line.FILTER.lower():
                print("Failing line for filter %s" % filter)
                continue
        if(line.REF == "C"):
            CTFreqArray.append(float(
                line.InfoDict["MAFS"].split(",")[2].split(">")[1]))
        elif(line.REF == "G"):
            GAFreqArray.append(
                line.InfoDict["MAFS"].split(",")[0].split(">")[1])
    pl("Length of C->T array: %s" % len(CTFreqArray), level=logging.DEBUG)
    pl("Length of G->A array: %s" % len(GAFreqArray), level=logging.DEBUG)
    GAFreqNP = np.array(GAFreqArray, dtype=np.longdouble)
    GAFreqNP = GAFreqNP[GAFreqNP < maxFreq]
    CTFreqNP = np.array(CTFreqArray, dtype=np.longdouble)
    CTFreqNP = CTFreqNP[CTFreqNP < maxFreq]
    GAFreqNP = GAFreqNP.reshape(-1, 1)
    CTFreqNP = CTFreqNP.reshape(-1, 1)
    if(concatenate):
        return np.concatenate([GAFreqNP, CTFreqNP])
    return GAFreqNP, CTFreqNP

