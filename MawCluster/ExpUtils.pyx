# cython: c_string_type=str, c_string_encoding=ascii
# cython: profile=True, cdivision=True, cdivision_warnings=True

import operator
import logging

import cython
import numpy as np
cimport numpy as np
from sklearn.covariance import EllipticEnvelope
from sklearn import svm
from scipy.misc import comb

from MawCluster.BCVCF import IterativeVCFFile
from utilBMF.HTSUtils import printlog as pl
from utilBMF.ErrorHandling import IllegalArgumentError
from MawCluster.SNVUtils import HeaderFilterDict

ctypedef np.float64_t dtype_t

"""
Contains tools for accurately identifying somatic mutations based on
different experimental protocols.
"""


# PROTOCOLS is a list of lowered strings.
PROTOCOLS = ["ffpe", "amplicon", "cf", "other"]


@cython.locals(FSD=cython.float, maxFreq=cython.float,
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
    cdef np.ndarray[dtype_t, ndim = 1] GAFreqNP
    cdef np.ndarray[dtype_t, ndim = 1] CTFreqNP
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
    GAFreqNP = np.array(GAFreqArray, dtype=np.float64)
    GAFreqNP = GAFreqNP[GAFreqNP < maxFreq]
    CTFreqNP = np.array(CTFreqArray, dtype=np.float64)
    CTFreqNP = CTFreqNP[CTFreqNP < maxFreq]
    GAFreqNP = GAFreqNP.reshape(-1, 1)
    CTFreqNP = CTFreqNP.reshape(-1, 1)
    if(concatenate):
        return np.concatenate([GAFreqNP, CTFreqNP])
    return GAFreqNP, CTFreqNP


@cython.locals(outliers_fraction=cython.float, contamination=cython.float,
               window=cython.long)
def DetectOutliers(f1, f2, outliers_fraction=0.1, contamination=0.005,
                   window=20):
    cdef np.ndarray[dtype_t, ndim = 1] GAFreqNP = f1
    cdef np.ndarray[dtype_t, ndim = 1] CTFreqNP = f2
    cdef np.ndarray[dtype_t, ndim = 1] FreqArray = np.concatenate(GAFreqNP,
                                                                  CTFreqNP)
    ee1 = EllipticEnvelope(contamination=contamination, assume_centered=False)
    ee2 = EllipticEnvelope(contamination=contamination, assume_centered=False)
    GAClassifier = ee1.fit(GAFreqNP)
    CTClassifier = ee2.fit(CTFreqNP)
    pass


def PlotNormalFit(array, outfile="default", maxFreq=0.2):
    import matplotlib as mpl 
    mpl.use('Agg')      # With this line = figure disappears
    import matplotlib.pyplot as plt 
    import matplotlib.mlab as mlab
    array = array.reshape(-1, 1)
    mu, sigma = np.mean(array), np.std(array)
    n, bins, patches = plt.hist(array, 50, normed=1, facecolor='green',
                                alpha=0.75)
    y = mlab.normpdf( bins, mu, sigma)
    l = plt.plot(bins, y, 'r--', linewidth=1)
    plt.xlabel('A->C/G->/T frequency')
    plt.ylabel('Proportion of sites')
    plt.title(r'$\mathrm{Histogram\ of\ deamination substitutions:}\ \mu=%s,\ \sigma=%s$' %(mu, sigma))
    plt.axis([np.min(array), np.max(array), 0., np.max(n)])
    plt.grid(True)
    plt.savefig(outfile + ".png")
    return outfile


@cython.locals(p=cython.float, k=cython.long, n=cython.long)
@cython.returns(cython.float)
def SamplingFrac(n, p=0., k=1):
    """
    Given a fixed probability of an event with n samplings, returns
    the probability that precisely "K" events have occurred.
    """
    assert 0. < p < 1. 
    return operator.mul(operator.pow(operator.sub(1, p), operator.sub(n, k)),
                        operator.pow(p, k)) * comb(n + 1, k)


@cython.locals(p=cython.float, k=cython.long, n=cython.long)
@cython.returns(cython.float)
def GetUnscaledProbs(n, p=0.):
    cdef np.ndarray[dtype_t, ndim = 1] arr
    arr = np.array([SamplingFrac(n, p=p, k=k) for k in range(n + 1)],
                    dtype=np.float64)
    return arr


@cython.locals(p=cython.float, k=cython.long, n=cython.long)
@cython.returns(cython.float)
def PartitionFunction(n, p=0.1):
    return np.sum(GetUnscaledProbs(n, p=p), 0)


@cython.locals(p=cython.float, n=cython.long, k=cython.long)
@cython.returns(cython.float)
def SamplingProbDist(n, p=0.):
    """
    Given a fixed probability of an event with n samplings, returns
    the probability that precisely "K" events have occurred.
    """
    assert 0. < p < 1. 
    cdef np.ndarray[dtype_t, ndim = 1] ProbDist
    cdef np.ndarray[dtype_t, ndim = 1] UnscaledProbs
    UnscaledProbs = GetUnscaledProbs(n, p=p)
    SummedUnscaledProbs = np.sum(UnscaledProbs, 0)
    ProbDist = np.divide(UnscaledProbs, SummedUnscaledProbs)
    return ProbDist


@cython.locals(maxPValue=cython.float)
def EstablishDeaminationBounds(inVCF, maxPValue=0.001, FILTER="default",
                               maxFreq=0.05):
    cdef np.ndarray[dtype_t, ndim = 1] FreqArray
    if(FILTER != "default"):
        FreqArray = GetDeaminationFrequencies(inVCF, maxFreq=maxFreq, FILTER=FILTER)
    else:
        FreqArray = GetDeaminationFrequencies(inVCF, maxFreq=maxFreq)
    AssumedProb = np.mean(FreqArray)
    


@cython.locals(numSD=cython.long)
def FilterByDeaminationFreq(inVCF, numSD=3, freqArray="default"):
    """
    Modifies the FILTER field of a VCF and adds a line to the header
    in an attempt to remove cytosine deamination artefacts.
    
    If the C->T or G->A frequency is more than `numSD` standard deviations
    above the mean, then that variant is not filtered out.
    Otherwise, DeaminationNoise replaces PASS or is appended to other filters.
    """
    
    
# somehow I need to work this into the analysis.
# apply the primer filter to amplicon