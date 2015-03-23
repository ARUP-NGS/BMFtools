# cython: c_string_type=str, c_string_encoding=ascii
# cython: profile=True, cdivision=True, cdivision_warnings=True

import cython
import numpy as np
cimport numpy as np
from scipy.stats import binom

from utilBMF.HTSUtils import printlog as pl, ThisIsMadness

ctypedef np.longdouble_t dtype128_t

"""
Contains probabilistic tools for accurately identifying somatic mutations.
"""

defaultPValue = 0.001

# PROTOCOLS is a list of lowered strings.
PROTOCOLS = ["ffpe", "amplicon", "cf", "other"]

@cython.locals(DOC=cython.long, pVal=dtype128_t,
               AC=cython.long)
def ConfidenceIntervalAAF(AC, DOC, pVal=defaultPValue):
    """
    Returns the confidence interval for an AAF given an observed AC
    and DOC.
    """
    return (np.array(binom.interval(1 - pVal, DOC,
                                    AC * 1. / DOC)) / DOC).tolist()


@cython.locals(n=cython.long, p=dtype128_t,
               pVal=dtype128_t)
@cython.returns(cython.long)
def GetCeiling(n, p=0.0, pVal=defaultPValue):
    """
    Returns the maximum fraction of events per sample with a p value of pVal,
    n samplings, and an assumed probability of p per sampling.
    """
    return binom.interval(1 - pVal, n, p)[1]
