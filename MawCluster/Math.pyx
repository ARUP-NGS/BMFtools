# cython: boundscheck=False, wraparound=False, initializedcheck=False
"""
Gamma functions for calculating Fisher scores.
"""
import cython

DEF LOG10E_X5_INV = 0.46051701859880917
# Multiply a phred score by this to convert

cdef inline double_t CHI2_FROM_PHRED(int32_t phredInt) nogil:
    return phredInt * LOG10E_X5_INV


cdef inline double_t INV_CHI2_FROM_PHRED(int32_t phredInt) nogil:
    return -2 * c_log(1 - c_pow(10, (phredInt * -1.)))


cdef inline double igamc_pvalues(int num_pvalues, double x) nogil:
    return 1.0 if(x < 0) else igamc(num_pvalues * 1., x / 2.0)