import cython
cimport cython
cimport numpy as np
from libc.stdint cimport int8_t, int32_t, int64_t, int16_t, uint64_t
from libc.stdio cimport fprintf, stderr as c_stderr
from utilBMF.Inliners cimport CONTEXT_TO_ARRAY_POS

cimport pysam.calignedsegment
cimport pysam.calignmentfile
ctypedef pysam.calignedsegment.AlignedSegment AlignedSegment_t
ctypedef pysam.calignmentfile.AlignmentFile AlignmentFile_t
ctypedef cython.str cystr
ctypedef double double_t
cimport cpython.array as c_array
ctypedef c_array.array py_array
from numpy cimport ndarray


@cython.boundscheck(False)
@cython.initializedcheck(False)
@cython.wraparound(False)
cdef inline int8_t NUC_TO_ARRAY(int8_t index) nogil:
    if(index == 0):
        return 65
    elif(index == 1):
        return 67
    elif(index == 2):
        return 71
    elif(index == 3):
        return 84
    else:
        return 78
