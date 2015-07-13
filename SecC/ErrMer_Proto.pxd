cimport cython
cimport pysam.calignedsegment
ctypedef pysam.calignedsegment.AlignedSegment AlignedSegment_t
ctypedef cython.str cystr

cdef inline cystr GetBS(AlignedSegment_t read):
    return read.opt("BS")
