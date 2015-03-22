# cython: c_string_type=str, c_string_encoding=ascii
# cython: profile=True, cdivision=True, cdivision_warnings=True


import logging
import operator
from operator import add as oadd
from operator import sub as osub
from operator import div as odiv
from operator import mul as omul
import math
from math import pow as mpow
from math import sqrt as msqrt

import cython
import numpy as np
cimport numpy as np
from scipy.misc import comb

from utilBMF.HTSUtils import printlog as pl, ThisIsMadness

ctypedef np.longdouble_t dtype128_t

"""
Contains probabilistic tools for accurately identifying somatic mutations.
"""

defaultPValue = 0.001


# PROTOCOLS is a list of lowered strings.
PROTOCOLS = ["ffpe", "amplicon", "cf", "other"]


@cython.locals(n=dtype128_t)
@cython.returns(dtype128_t)
def StirlingsApprox(n):
    """
    Stirling's Approximation is a continuous function which approximates
    factorials extremely accurately.
    n! ~ Sqrt(2*Pi*n) * (n / e)^n
    """
    return omul(mpow(odiv(n, np.e), n), msqrt(reduce(omul, [2, np.pi, n], 1)))


@cython.locals(n=dtype128_t, k=dtype128_t)
@cython.returns(dtype128_t)
def StirlingsFact(n, k):
    """
    Stirling's Approximation is a continuous function which approximates
    factorials extremely accurately.
    n! ~ Sqrt(2*Pi*n) * (n / e)^n
    """
    return StirlingsApprox(n) / StirlingsApprox(k) / StirlingsApprox(
        osub(n, k))


@cython.locals(p=dtype128_t, k=cython.long, n=cython.long)
@cython.returns(dtype128_t)
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


@cython.locals(p=dtype128_t, k=cython.long, n=cython.long)
@cython.returns(dtype128_t)
def SamplingFrac_(n, p=0., k=1):
    """
    Given a fixed probability of an event with n samplings, returns
    the probability moment (IE, unnormalized) that precisely "k" events
    Uses Stirling's Approximation for speed.
    """
    assert 0. < p < 1.
    return reduce(omul, [mpow(osub(1, p), osub(n, k)), mpow(p, k),
                         StirlingsFact(n + 1, k)], 1)


@cython.locals(p=dtype128_t, k=cython.long, n=cython.long)
@cython.returns(np.ndarray)
def GetUnscaledProbs(n, p=0.):
    """
    Calculates the probability moments for k in range(n + 1).
    This uses a factorial, primarily for
    speed, but it's also definitely good enough.
    """
    return np.array([SamplingFrac(n, p=p, k=k) for k in range(n + 1)],
                    dtype=np.longdouble)


@cython.locals(p=dtype128_t, k=cython.long, n=cython.long)
@cython.returns(np.ndarray)
def GetUnscaledProbs_(n, p=0.):
    """
    Calculates the probability moments for k in range(n + 1).
    This uses Stirling's approximation instead of factorial, primarily for
    speed, but it's also definitely good enough.
    """
    return np.array([SamplingFrac_(n, p=p, k=k) for k in range(n + 1)],
                    dtype=np.longdouble)


@cython.locals(p=dtype128_t, k=cython.long, n=cython.long)
@cython.returns(dtype128_t)
def PartitionFunction(n, p=0.1):
    return np.sum(GetUnscaledProbs(n, p=p), 0)


@cython.locals(p=dtype128_t, k=cython.long, n=cython.long)
@cython.returns(dtype128_t)
def PartitionFunction_(n, p=0.1):
    return np.sum(GetUnscaledProbs_(n, p=p), 0)


@cython.locals(p=dtype128_t, n=cython.long, k=cython.long)
@cython.returns(np.ndarray)
def SamplingProbDist(n, p=0.):
    """
    Given a fixed probability of an event with n samplings, returns
    the probability that precisely "K" events have occurred.
    """
    assert 0. < p < 1.
    cdef np.ndarray[dtype128_t, ndim = 1] ProbDist
    cdef np.ndarray[dtype128_t, ndim = 1] UnscaledProbs
    UnscaledProbs = GetUnscaledProbs(n, p=p)
    PartitionFn = np.sum(UnscaledProbs, 0)
    ProbDist = UnscaledProbs / PartitionFn
    return ProbDist


@cython.locals(p=dtype128_t, n=cython.long, k=cython.long)
@cython.returns(np.ndarray)
def SamplingProbDist_(n, p=0.):
    """
    Given a fixed probability of an event with n samplings, returns
    the probability that precisely "K" events have occurred.
    """
    assert 0. < p < 1.
    cdef np.ndarray[dtype128_t, ndim = 1] ProbDist
    cdef np.ndarray[dtype128_t, ndim = 1] UnscaledProbs
    UnscaledProbs = GetUnscaledProbs_(n, p=p)
    PartitionFn = np.sum(UnscaledProbs, 0)
    ProbDist = UnscaledProbs / PartitionFn
    return ProbDist


@cython.locals(p=dtype128_t, n=cython.long, k=cython.long,
               PartitionFn=dtype128_t)
def SamplingProbMoments(n, p=0.):
    """
    Given a fixed probability of an event with n samplings, returns
    the probability that precisely "K" events have occurred.
    """
    assert 0. < p < 1.
    cdef np.ndarray[dtype128_t, ndim = 1] UnscaledProbs
    UnscaledProbs = GetUnscaledProbs(n, p=p)
    PartitionFn = np.sum(UnscaledProbs, 0, dtype=np.longdouble)
    return UnscaledProbs, PartitionFn


@cython.locals(p=dtype128_t, n=cython.long, k=cython.long,
               PartitionFn=dtype128_t)
def SamplingProbMoments_(n, p=0.):
    """
    Given a fixed probability of an event with n samplings, returns
    the probability that precisely "K" events have occurred.
    """
    assert 0. < p < 1.
    cdef np.ndarray[dtype128_t, ndim = 1] UnscaledProbs
    UnscaledProbs = GetUnscaledProbs_(n, p=p)
    PartitionFn = np.sum(UnscaledProbs, 0, dtype=np.longdouble)
    return UnscaledProbs, PartitionFn


@cython.locals(p=dtype128_t, pVal=dtype128_t,
               n=cython.long, maxProb=cython.long, minFreq=dtype128_t,
               maxFreq=dtype128_t, scaledPVal=dtype128_t,
               PartitionFn=dtype128_t)
def ConfidenceInterval(n, p=0.0, pVal=defaultPValue):
    """
    Returns the minimum and maximum frequency for an event taking place at
    probability p, the number of samplings n, and a pVal of (1-conf), where
    conf is the fraction confidence interval desired.
    """
    if(pVal == defaultPValue):
        pl("pVal has not been changed. Default of 0.1% is being used.",
           level=logging.DEBUG)
    cdef np.ndarray[dtype128_t, ndim = 1] ProbMoments
    ProbMoments, PartitionFn = SamplingProbMoments(n, p=p)
    scaledPVal = pVal * PartitionFn
    maxProb = np.argmax(ProbMoments)
    for i in range(maxProb):
        if(operator.ge(oadd(sum(
                ProbMoments[osub(maxProb, i):oadd(maxProb, i)]),
                                    scaledPVal), PartitionFn)):
            minFreq = odiv(osub(maxProb, i), float(n))
            maxFreq = odiv(oadd(maxProb, i), float(n))
            return minFreq, maxFreq, pVal
    minFreq = -1.
    maxFreq = -1.
    pl("Confidence interval could not meet the desired p value.",
       level=logging.DEBUG)
    return minFreq, maxFreq, pVal


@cython.locals(p=dtype128_t, pVal=dtype128_t,
               n=cython.long, maxProb=cython.long, minFreq=dtype128_t,
               maxFreq=dtype128_t, scaledPVal=dtype128_t,
               PartitionFn=dtype128_t)
def ConfidenceInterval_(n, p=0.0, pVal=defaultPValue):
    """
    Returns the minimum and maximum frequency for an event taking place at
    probability p, the number of samplings n, and a pVal of (1-conf), where
    conf is the fraction confidence interval desired.
    """
    if(pVal == defaultPValue):
        pl("pVal has not been changed. Default of 0.1% is being used.",
           level=logging.DEBUG)
    cdef np.ndarray[dtype128_t, ndim = 1] ProbMoments
    ProbMoments, PartitionFn = SamplingProbMoments_(n, p=p)
    scaledPVal = pVal * PartitionFn
    maxProb = np.argmax(ProbMoments)
    for i in range(maxProb):
        if(operator.ge(oadd(sum(
                ProbMoments[osub(maxProb, i):oadd(maxProb, i)]),
                                    scaledPVal), PartitionFn)):
            minFreq = odiv(osub(maxProb, i), float(n))
            maxFreq = odiv(oadd(maxProb, i), float(n))
            return minFreq, maxFreq, pVal
    minFreq = -1.
    maxFreq = -1.
    pl("Confidence interval could not meet the desired p value.",
       level=logging.DEBUG)
    return minFreq, maxFreq, pVal

