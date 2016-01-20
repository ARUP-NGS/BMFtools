from libc.stdlib cimport malloc, free
from libc.string cimport strcmp, memcpy
from libc.stdint cimport *
from cpython.string cimport PyString_AsString
from libc.stdlib cimport EXIT_FAILURE, qsort
cimport cython
from cython cimport view
from numpy cimport ndarray, npy_intp, int64_t, uint8_t
from cpython cimport array as c_array
ctypedef c_array.array py_array
ctypedef cython.str cystr
cimport numpy as np

cpdef py_array str2intarray(cystr instr)

cdef inline py_array cs_to_ia(cystr input_str)
cdef inline py_array cs_to_ph(cystr input_str)


cdef public cystr PH2CHR_TRANS

cdef char *revcmp(char *, char *, uint64_t l)

cdef struct struct_str:
    char * string
    size_t size
