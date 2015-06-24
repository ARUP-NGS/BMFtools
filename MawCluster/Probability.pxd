cimport cython
cimport numpy as np
from numpy cimport ndarray
from cpython cimport array as c_array
ctypedef c_array.array py_array
ctypedef long double npy_float128
ctypedef npy_float128   float128_t
ctypedef double npy_float64
ctypedef npy_float64   float64_t
ctypedef cython.str cystr

from libc.math cimport sqrt as c_sqrt, M_SQRT1_2, log as c_log

cdef inline double square(double input) nogil:
    return input * input

cdef inline double c_abs(double input) nogil:
    if(input < 0):
        return -1 * input
    else:
        return input

cpdef float64_t HellingerDistance(ndarray[float64_t, ndim=1] array1,
                                  ndarray[float64_t, ndim=1] array2)
cdef float64_t cHellingerDistance(float64_t* array1,
                                  float64_t* array2,
                                  size_t length)
cpdef float64_t HellingDistanceDict(dict dict1, dict dict2)