from libc.stdlib cimport malloc, free
from libc.string cimport strcmp
from cpython.string cimport PyString_AsString
cimport cython
from cython cimport view
from numpy cimport ndarray, npy_intp, int64_t, uint8_t
from cpython cimport array as c_array
from utilBMF.Inliners cimport RevCmpChar, RevCmpInt, RevCmpToChar
ctypedef c_array.array py_array
ctypedef cython.str cystr
cimport numpy as np

cdef py_array cs_to_ia(cystr input_str)
cpdef py_array str2intarray(cystr instr)
cdef py_array cs_to_ph(cystr input_str)

cdef cystr RevCmpImplicit(cystr seq)

