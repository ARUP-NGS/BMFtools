# Math functions

from libc.math cimport exp, INFINITY, M_PI as PI
from libc.math cimport floor, sin, log, fabs, log as c_log, pow as c_pow
from libc.stdlib cimport malloc, free, realloc
from libc.stdint cimport int32_t
ctypedef double double_t


# FUNCTIONS
cdef inline double_t igamc(double_t a, double_t x) nogil
cdef inline double_t CHI2_FROM_PHRED(int32_t phredInt) nogil
cdef inline double_t INV_CHI2_FROM_PHRED(int32_t phredInt) nogil

DEF LOG10E_X5_INV = 0.46051701859880917
# Multiply a phred score by this to convert
