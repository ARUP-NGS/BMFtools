from numpy cimport ndarray
cimport numpy as np
from libc.stdlib cimport malloc, free
cdef extern from 'FqProc.h':
    ctypedef struct Read:
        np.int8_t* sequence
        np.int8_t* quality
        np.int8_t length;
    cdef Read FqMerge (int** sequence, int** quality, int nreads, int lenreads)