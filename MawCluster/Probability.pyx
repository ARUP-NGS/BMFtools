# cython: c_string_type=str, c_string_encoding=ascii

import logging
import operator
from operator import add as oadd
from operator import sub as osub
from operator import div as odiv
from operator import mul as omul
import math
from math import pow as mpow
from math import sqrt as msqrt
from utilBMF.ErrorHandling import ThisIsMadness

import cython
import numpy as np
from numpy import mean as nmean
from numpy import sum as nsum
from numpy cimport ndarray
from scipy.stats import binom
from scipy.misc import comb
from statsmodels.stats import proportion
from cytoolz import memoize
pconfint = proportion.proportion_confint



"""
Contains probabilistic tools for accurately identifying somatic mutations.
"""

defaultPValue = 0.05

# PROTOCOLS is a list of lowered strings.
PROTOCOLS = ["ffpe", "amplicon", "cf", "other"]


cpdef float64_t HellingDistanceDict(dict dict1, dict dict2):
    cdef cystr key
    cdef ndarray[float64_t, ndim=1] array1, array2
    array1 = np.array([dict1[key] for key in
                       sorted(dict1.iterkeys())])
    array2 = np.array([dict2[key] for key in
                       sorted(dict2.iterkeys())])
    return HellingerDistance(array1, array2)


cpdef float64_t HellingerDistance(ndarray[float64_t, ndim=1] array1,
                                  ndarray[float64_t, ndim=1] array2):
    cdef size_t length = len(array1)
    assert length == len(array2)  # Sanity check
    return cHellingerDistance(&array1[0], &array2[0], length)


cpdef float64_t BhattacharyyaDistance(ndarray[float64_t, ndim=1] array1,
                                      ndarray[float64_t, ndim=1] array2):
    cdef size_t length = len(array1)
    assert length == len(array2)  # Sanity check
    return cBhattacharyyaDistance(&array1[0], &array2[0], length)


cdef float64_t cBhattacharyyaDistance(float64_t* array1,
                                      float64_t* array2,
                                      size_t length):
    """
    Calculates the Helling Distance between two discrete probability
    distributions, described (each) as a 1-dimensional array.
    """
    cdef float64_t cumSum = 0.
    cdef float64_t tmpFloat
    cdef size_t tmpInt
    for tmpInt in range(length):
        cumSum += c_sqrt(array1[tmpInt] * array2[tmpInt])
    return -1 * c_log(cumSum)


cdef float64_t cHellingerDistance(float64_t* array1,
                                  float64_t* array2,
                                  size_t length):
    """
    Calculates the Helling Distance between two discrete probability
    distributions, described (each) as a 1-dimensional array.
    """
    cdef float64_t cumSum = 0.
    cdef float64_t tmpFloat
    cdef size_t tmpInt
    for tmpInt in range(length):
        cumSum += c_abs(c_sqrt(array1[tmpInt]) -
                        c_sqrt(array2[tmpInt]))
    return cumSum * M_SQRT1_2


@memoize
@cython.locals(DOC=int, pVal=float128_t,
               AC=int)
def ConfidenceIntervalAAF(AC, DOC, pVal=defaultPValue,
                          method="agresti_coull"):
    """
    Returns the confidence interval for an AAF given an observed AC
    and DOC.
    """
    if(method == "scipy"):
        return (np.array(binom.interval(1 - pVal, DOC,
                                       AC * 1. / DOC)) / DOC).tolist()
    try:
        return pconfint(AC, DOC, alpha=pVal, method=method)
    except NotImplementedError:
        raise ThisIsMadness("Confidence interval method `%s` not " % method +
                            "implemented! Check the docs for statsmodels."
                            "stats.proportion.")


@cython.returns(ndarray)
def ConfidenceIntervalAI(int Allele1,
                         int Allele2,
                         float128_t pVal=defaultPValue,
                         cystr method="agresti_coull"):
    """
    Returns the confidence interval for an allelic imbalance
    given counts for Allele1 and Allele2, where those are the most common
    alleles at a given position.
    """
    cdef float128_t ratio
    if(Allele1 < Allele2):
        ratio = 1. * Allele1 / Allele2
    else:
        ratio = 1. * Allele2 / Allele1
    if(method == "scipy"):
        return np.array(binom.interval(1 - pVal,
                                      Allele1 + Allele2, ratio),
                       dtype=np.float128)
    try:
        return np.array(pconfint(Allele1, Allele1 + Allele2,
                                alpha=pVal, method=method),
                       dtype=np.float128)
    except NotImplementedError:
        raise ThisIsMadness("Confidence interval method `%s` not" % method +
                            " implemented! Check the docs for statsmodels."
                            "stats.proportion.")


@cython.returns(tuple)
def MakeAICall(int Allele1,
               int Allele2,
               float128_t pVal=defaultPValue,
               cystr method="agresti_coull"):
    """
    Gets confidence bounds, returns a call, a "minor" allele frequency,
    an observed allelic imbalance ratio (as defined by more common allele
    over less common allele), and a confidence interval
    """
    cdef float128_t minorAF
    cdef cython.bint call
    cdef float128_t allelicImbalanceRatio
    confInt = ConfidenceIntervalAI(Allele1, Allele2, pVal=defaultPValue,
                                   method=method)
    minorAF = nmean(confInt, dtype=np.float128)
    allelicImbalanceRatio = 1. * minorAF / (1 - minorAF)
    call = (minorAF < confInt[0] or minorAF > confInt[1])
    return call, allelicImbalanceRatio, confInt


@memoize
@cython.locals(n=int, p=float128_t,
               pVal=float128_t)
@cython.returns(int)
def GetCeiling(n, p=0.0, pVal=defaultPValue):
    """
    Returns the maximum fraction of events per sample with a p value of pVal,
    n samplings, and an assumed probability of p per sampling.
    """
    return binom.interval(1 - pVal, n, p)[1]


@memoize
@cython.locals(n=float128_t)
@cython.returns(float128_t)
def StirlingsApprox(n):
    """
    Stirling's Approximation is a continuous function which approximates
    factorials extremely accurately.
    n! ~ Sqrt(2*Pi*n) * (n / e)^n
    """
    return omul(mpow(odiv(n, np.e), n), msqrt(reduce(omul, [2, np.pi, n], 1)))


@memoize
@cython.locals(n=float128_t, k=float128_t)
@cython.returns(float128_t)
def StirlingsFact(n, k):
    """
    Stirling's Approximation is a continuous function which approximates
    factorials extremely accurately.
    n! ~ Sqrt(2*Pi*n) * (n / e)^n
    """
    return StirlingsApprox(n) / StirlingsApprox(k) / StirlingsApprox(
        osub(n, k))


@memoize
@cython.locals(p=float128_t, k=int, n=int)
@cython.returns(float128_t)
def SamplingFrac(n, p=0., k=1):
    """
    Given a fixed probability of an event with n samplings, returns
    the probability moment (IE, unnormalized) that precisely "k" events
    have occurred.
    If you want speed, use SamplingFrac_, which uses Stirling's Approximation.
    """
    assert 0. < p < 1.
    return reduce(omul, [mpow(osub(1, p), osub(n, k)),
                  mpow(p, k), comb(n + 1, k)], 1)


@memoize
@cython.locals(p=float128_t, k=int, n=int)
@cython.returns(float128_t)
def SamplingFrac_(n, p=0., k=1):
    """
    Given a fixed probability of an event with n samplings, returns
    the probability moment (IE, unnormalized) that precisely "k" events
    Uses Stirling's Approximation for speed.
    """
    assert 0. < p < 1.
    return reduce(omul, [mpow(osub(1, p), osub(n, k)), mpow(p, k),
                  StirlingsFact(n + 1, k)], 1)


@memoize
@cython.locals(p=float128_t, k=int, n=int)
@cython.returns(ndarray)
def GetUnscaledProbs(n, p=0.):
    """
    Calculates the probability moments for k in xrange(n + 1).
    This uses a factorial, primarily for
    speed, but it's also definitely good enough.
    """
    return np.array([SamplingFrac(n, p=p, k=k) for k in xrange(n + 1)],
                   dtype=np.float128)


@memoize
@cython.locals(p=float128_t, k=int, n=int)
@cython.returns(ndarray)
def GetUnscaledProbs_(n, p=0.):
    """
    Calculates the probability moments for k in xrange(n + 1).
    This uses Stirling's approximation instead of factorial, primarily for
    speed, but it's also definitely good enough.
    """
    return np.array([SamplingFrac_(n, p=p, k=k) for k in xrange(n + 1)],
                   dtype=np.float128)


@memoize
@cython.locals(p=float128_t, k=int, n=int)
@cython.returns(float128_t)
def PartitionFunction(n, p=0.1):
    return nsum(GetUnscaledProbs(n, p=p), 0)


@memoize
@cython.locals(p=float128_t, k=int, n=int)
@cython.returns(float128_t)
def PartitionFunction_(n, p=0.1):
    return nsum(GetUnscaledProbs_(n, p=p), 0)


@memoize
@cython.locals(p=float128_t, n=int, k=int)
@cython.returns(ndarray)
def SamplingProbDist(n, p=0.):
    """
    Given a fixed probability of an event with n samplings, returns
    the probability that precisely "K" events have occurred.
    """
    assert 0. < p < 1.
    cdef ndarray[float128_t, ndim=1] ProbDist
    cdef ndarray[float128_t, ndim=1] UnscaledProbs
    UnscaledProbs = GetUnscaledProbs(n, p=p)
    PartitionFn = nsum(UnscaledProbs, 0)
    ProbDist = UnscaledProbs / PartitionFn
    return ProbDist


@memoize
@cython.locals(p=float128_t, n=int, k=int)
@cython.returns(ndarray)
def SamplingProbDist_(n, p=0.):
    """
    Given a fixed probability of an event with n samplings, returns
    the probability that precisely "K" events have occurred.
    """
    assert 0. < p < 1.
    cdef ndarray[float128_t, ndim=1] ProbDist
    cdef ndarray[float128_t, ndim=1] UnscaledProbs
    UnscaledProbs = GetUnscaledProbs_(n, p=p)
    PartitionFn = nsum(UnscaledProbs, 0)
    ProbDist = UnscaledProbs / PartitionFn
    return ProbDist


@memoize
@cython.locals(p=float128_t, n=int, k=int,
               PartitionFn=float128_t)
def SamplingProbMoments(n, p=0.):
    """
    Given a fixed probability of an event with n samplings, returns
    the probability that precisely "K" events have occurred.
    """
    assert 0. < p < 1.
    cdef ndarray[float128_t, ndim=1] UnscaledProbs
    UnscaledProbs = GetUnscaledProbs(n, p=p)
    PartitionFn = nsum(UnscaledProbs, 0, dtype=np.float128)
    return UnscaledProbs, PartitionFn


@memoize
@cython.locals(p=float128_t, n=int, k=int,
               PartitionFn=float128_t)
def SamplingProbMoments_(n, p=0.):
    """
    Given a fixed probability of an event with n samplings, returns
    the probability that precisely "K" events have occurred.
    """
    assert 0. < p < 1.
    cdef ndarray[float128_t, ndim=1] UnscaledProbs
    UnscaledProbs = GetUnscaledProbs_(n, p=p)
    PartitionFn = nsum(UnscaledProbs, 0, dtype=np.float128)
    return UnscaledProbs, PartitionFn
