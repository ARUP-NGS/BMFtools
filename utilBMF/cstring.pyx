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


cpdef print_chars1(cystr input_str):
    """Looks like it's 3x as fast to use the direct C conversion
    rather than making the views and accessing them, based on my comparison
    with print_chars
    """
    cdef carray arr = cs_to_ia(input_str)
    cdef int i
    for i in xrange(100000):
        ",".join(map(str, arr))
    return


cdef carray cs_to_ia2(cystr input_str):
    return array('i', ps2va(<char*> input_str, len(input_str)))


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef carray str2intarray2(cystr input_str):
    """Don't use this function - it is equivalent to str2intarray, but
    slower than naive python
    """
    return <carray> ps2va(input_str, len(input_str))


cdef carray cs_to_ia(cystr input_str):
    cdef char i
    return array('B', [i for i in <char *> input_str])


cdef carray cs_to_ph(cystr input_str):
    cdef char i
    return array('B', [i - 33 for i in <char *> input_str])


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef carray str2intarray(cystr instr):
    return cs_to_ia(instr)
