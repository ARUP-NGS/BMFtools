cimport cython
cimport pysam.calignmentfile
ctypedef pysam.calignmentfile.AlignedSegment AlignedSegment_t
ctypedef cython.str cystr

cdef inline cystr GetBS(AlignedSegment_t read):
    return read.opt("BS")
