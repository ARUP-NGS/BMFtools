cimport cython
cimport numpy as np
ctypedef double npy_float64
ctypedef npy_float64   float64_t
ctypedef cython.str cystr


cdef extern from "../dlib/hellinger.c":
    cdef float64_t Hellinger_in_c(
        float64_t* arr1, float64_t* arr2, size_t length)
