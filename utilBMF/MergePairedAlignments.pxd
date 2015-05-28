cimport cython
cimport pysam.calignmentfile

cdef class LayoutPos(object):
    cdef public cython.int pos, readPos, quality, agreement
    cdef public cython.str operation, base

cdef class Layout(object):
    cdef public pysam.calignmentfile.AlignedSegment read
    cdef public list positions
    cdef public cython.int firstMapped