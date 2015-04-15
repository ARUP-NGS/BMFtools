cimport cython
cimport pysam.calignmentfile

cdef class pPileupRead:
    """
    Python container for the PileupRead proxy in pysam
    """
    cdef public cython.str BaseCall
    cdef public cython.bint is_del
    cdef public cython.long level
    cdef public cython.long indel
    cdef public cython.long query_position
    cdef public cython.str name
    cdef public pysam.calignmentfile.AlignedSegment alignment
