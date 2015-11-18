from libc.stdlib cimport malloc, free
from libc.string cimport strcmp, memcpy
from cpython.string cimport PyString_AsString
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

cdef struct struct_str:
    char * string
    size_t size

cdef extern from "../lib/cstr_util.h" nogil:
    int lex_lt(char *s, size_t l)
    char *revcmp(char *, char *, size_t)