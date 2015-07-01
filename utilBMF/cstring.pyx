# @cython: c_string_encoding=ascii
import numpy as np
from numpy import frombuffer
import cython
import ctypes
from array import array


cdef int * to_cstring_array(cystr input_str):
    cdef char *tmp = PyString_AsString(input_str)
    cdef int *ret = <int *>malloc(len(input_str) * sizeof(char *))
    ret = <int*> tmp
    return ret

# define a function that can deallocate the data (if needed)
# my_array.callback_free_data = free


cdef view.array ps2va(char * inArr, size_t size):
    """This code is here until I figure out how to work with the string's values
    directly as integers without having to copy. This somehow seems close but
    just isn't working as I'd want it to.
    """
    cdef view.array ret_array = view.array(shape=(size,),
                                           itemsize=sizeof(char),
                                           format="i", allocate_buffer=False,
                                           mode='c')
    ret_array.data = inArr
    return ret_array


cpdef print_chars(cystr input_str):
    """Looks like it's 3x as fast to use the direct C conversion
    rather than making the views and accessing them, based on my comparison
    with print_chars
    """
    cdef py_array arr = cs_to_ia(input_str)
    cdef int i
    for i in xrange(100000):
        ",".join(map(str, arr))
    return


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.initializedcheck(False)
cdef py_array cs_to_ia(cystr input_str):
    cdef char i
    return array('B', input_str)


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.initializedcheck(False)
cdef inline cystr RevCmpPyArray(cystr seq):
    cdef char i
    return array('B', [RevCmpToChar(i) for
                       i in <char *>seq]).tostring()[::-1]


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.initializedcheck(False)
cdef cystr RevCmpImplicit(cystr seq):
    """
    Very quickly reverse complements a string with an inline switch for faster
    memoization than dictionary access.
    """
    cdef char i
    return "".join([RevCmpInt(i) for i in <char *>seq])[::-1]


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef py_array str2intarray(cystr instr):
    return cs_to_ia(instr)


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef py_array str2phredarray(cystr instr):
    return cs_to_ph(instr)
