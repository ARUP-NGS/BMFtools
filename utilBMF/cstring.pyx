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


cdef py_array cs_to_ia2(cystr input_str):
    return array('i', ps2va(<char*> input_str, len(input_str)))


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.initializedcheck(False)
cpdef py_array str2intarray2(cystr input_str):
    """Don't use this function - it is equivalent to str2intarray, but
    slower than naive python
    """
    return <py_array> ps2va(input_str, len(input_str))


cdef py_array cs_to_ia(cystr input_str):
    cdef char i
    return array('B', [i for i in <char *> input_str])


'''
Speed experiments

In [5]: %timeit strlen0("Hello")
The slowest run took 46.08 times longer than the fastest. This could mean that an intermediate result is being cached
10000000 loops, best of 3: 88 ns per loop

In [6]: %timeit strlen2("Hello")
The slowest run took 68.72 times longer than the fastest. This could mean that an intermediate result is being cached
10000000 loops, best of 3: 59 ns per loop

In [7]: %timeit strlen1("Hello")
The slowest run took 63.71 times longer than the fastest. This could mean that an intermediate result is being cached
10000000 loops, best of 3: 48.7 ns per loop

In [8]: %timeit len("Hello")
The slowest run took 74.42 times longer than the fastest. This could mean that an intermediate result is being cached
10000000 loops, best of 3: 41.6 ns per loop

It makes sense that strlen1 would be slower than
'''

cdef int cstrlen0(char * input_str):
    cdef int count = 0
    cdef char i
    for i in input_str:
        count += 1
    return count

cpdef int strlen0(cystr input_str):
    return cstrlen0(<char *>input_str)

cdef int cstrlen1(char * input_str, int length):
    cdef int count = 0
    cdef char i
    for char in range(length):
        count += 1
    return count

cpdef int strlen1(cystr input_str):
    return cstrlen1(<char *>input_str, len(input_str))


cdef int cstrlen2(cystr input_str):
    cdef int count = 0
    cdef char i
    for char in <char*> input_str:
        count += 1
    return count

cpdef int strlen2(cystr input_str):
    return cstrlen2(input_str)
'''
Speed tests ---
In [26]: %timeit char_in_str("Hello", "o")
10000000 loops, best of 3: 79.9 ns per loop

In [27]: %timeit char_in_str1("Hello", "o")
10000000 loops, best of 3: 92.8 ns per loop

In [28]: %timeit char_in_cstr("Hello", "o")
10000000 loops, best of 3: 67.6 ns per loop

In [29]: %timeit char_in_p_str("Hello", "o")
10000000 loops, best of 3: 69.7 ns per loop

char_in_cstr is faster than char_in_p_str (naive python stdlib),
but only be about 3%, and even that uses the standard
python syntax of x in y.

'''

@cython.wraparound(False)
@cython.boundscheck(False)
@cython.nonecheck(False)
@cython.initializedcheck(False)
cdef bint char_in_cstr(bytes input_str, char char_arg):
    return char_arg in input_str

@cython.wraparound(False)
@cython.boundscheck(False)
@cython.nonecheck(False)
@cython.initializedcheck(False)
cpdef bint char_in_str2(cystr input_str, cystr char_arg):
    '''
    cdef char * char_arg_ptr = char_arg
    cdef char tmpInt = char_arg_ptr[0]
    '''
    return char_in_cstr(input_str,
                         (<char*> char_arg)[0])

@cython.wraparound(False)
@cython.boundscheck(False)
@cython.nonecheck(False)
@cython.initializedcheck(False)
cdef bint char_in_c_str1(char * input_str, char char_arg):
    return char_arg in input_str

@cython.wraparound(False)
@cython.boundscheck(False)
@cython.nonecheck(False)
@cython.initializedcheck(False)
cpdef bint char_in_str1(cystr input_str, cystr char_arg):
    '''
    cdef char * char_arg_ptr = char_arg
    cdef char tmpInt = char_arg_ptr[0]
    '''
    return char_in_c_str1(<char *>input_str,
                          (<char*> char_arg)[0])

@cython.wraparound(False)
@cython.boundscheck(False)
@cython.nonecheck(False)
@cython.initializedcheck(False)
cdef bint char_in_c_str_manual(char * input_str, int length, char char_arg):
    cdef char i
    for i in range(length):
        if(input_str[i] == char_arg):
            return True
    return False

@cython.wraparound(False)
@cython.boundscheck(False)
@cython.nonecheck(False)
@cython.initializedcheck(False)
cpdef bint char_in_str(cystr input_str, cystr char_arg):
    '''
    cdef char * char_arg_ptr = char_arg
    cdef char tmpInt = char_arg_ptr[0]
    '''
    return char_in_c_str_manual(<char *>input_str, len(input_str),
                                (<char*> char_arg)[0])

@cython.wraparound(False)
@cython.boundscheck(False)
@cython.nonecheck(False)
@cython.initializedcheck(False)
@cython.returns(bint)
@cython.locals(input_str=cystr, char_arg=cystr)
def char_in_p_str(input_str, char_arg):
    return char_arg[0] in input_str


cdef py_array cs_to_ph(cystr input_str):
    cdef char i
    return array('B', [i - 33 for i in <char *> input_str])


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef py_array str2intarray(cystr instr):
    return cs_to_ia(instr)


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef py_array str2phredarray(cystr instr):
    return cs_to_ph(instr)
