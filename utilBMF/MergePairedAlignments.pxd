cimport cython
cimport pysam.calignmentfile
ctypedef Layout Layout_t

cdef class LayoutPos(object):
    cdef public cython.int pos, readPos, quality, agreement
    cdef public cython.str operation, base
    cpdef cython.bint ismapped(self)

cdef class Layout(object):
    cdef public pysam.calignmentfile.AlignedSegment read
    cdef public list positions
    cdef public dict tagDict
    cdef public cython.int firstMapped, InitPos, flag, pnext, tlen, mapq
    cdef public cython.str Name, contig, rnext
    cdef public cython.bint isMerged, is_reverse

    cpdef cython.int getAlignmentStart(self)
    cpdef cython.str getCigarString(self)
    cpdef cython.str getSeq(self)

cpdef Layout_t MergeLayoutsToLayout(Layout_t, Layout_t)