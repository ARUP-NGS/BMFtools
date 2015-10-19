# Math functions

from libc.math cimport exp, INFINITY, M_PI as PI
from libc.math cimport floor, sin, log, fabs, log as c_log, pow as c_pow, log10 as c_log10
from libc.stdlib cimport malloc, free, realloc
from libc.stdint cimport int32_t, int8_t
ctypedef double double_t


# FUNCTIONS
cdef inline double_t CHI2_FROM_PHRED(int32_t phredInt) nogil
cdef inline double_t INV_CHI2_FROM_PHRED(int32_t phredInt) nogil
cdef int8_t * arrmax(int32_t * quals, int8_t * ret,
                     size_t nRecs, size_t rLen) nogil


cdef extern from "../include/igamc_cephes.h":
    cdef double igamc(double, double) nogil

cpdef inline double_t pigamc(double a, double b):
    return igamc(a, b)


cdef inline int MergeAgreedQualities(int q1, int q2) nogil:
    return <int>(-10 * c_log10(igamc(2., CHI2_FROM_PHRED(q1) +
                                     CHI2_FROM_PHRED(q2) / 2.0)) + 0.5)


cdef inline int MergeDiscQualities(int q1, int q2) nogil:
    if(q1 > q2):
        return <int>(- 10 * c_log10(igamc(2., INV_CHI2_FROM_PHRED(q2) +
                                          CHI2_FROM_PHRED(q1))) + 0.5)
    else:
        return <int>(- 10 * c_log10(igamc(2., INV_CHI2_FROM_PHRED(q1) +
                                          CHI2_FROM_PHRED(q2))) + 0.5)
