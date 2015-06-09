from libc.stdlib cimport malloc, free
from libc.string cimport strcmp, strlen
from cpython.string cimport PyString_AsString
cimport cython
from cython cimport view
from numpy cimport ndarray, npy_intp, int64_t, uint8_t
from cpython cimport array as c_array
ctypedef c_array.array carray
ctypedef cython.str cystr
cimport numpy as np

cdef carray cs_to_ia(cystr input_str)
cpdef carray str2intarray(cystr instr)
cdef carray cs_to_ph(cystr input_str)