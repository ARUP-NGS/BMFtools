# Math functions

ctypedef double double_t
from libc.math cimport exp, INFINITY, M_PI as PI, floor, sin, log, fabs
from libc.stdlib cimport malloc, free, realloc


# FUNCTIONS
cdef inline double_t igamc(double_t a, double_t x) nogil