cimport cython
cimport pysam.calignmentfile
cimport numpy as np
ctypedef Layout Layout_t
ctypedef LayoutPos LayoutPos_t

cdef class LayoutPos(object):
    cdef public cython.int pos, readPos, quality, agreement
    cdef public char operation, base
    # cdef public cython.str operation, base
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
    cdef np.ndarray[char] getSeqArr(self, dict chrDict=?)
    cdef int getRP4FP(self)
    cpdef int getRefPosForFirstPos(self)
    cpdef np.ndarray[int] getAgreement(self)
    cdef np.ndarray[int] getAgreementNP(self)
    cdef np.ndarray[int] getQualNP(self)
    cpdef np.ndarray[int] getQual(self)
    cdef cython.str getQualStringNP(self, dict ph2chrDict=?)
    cpdef getQualString(self)
    cdef int getLastRefPosNP(self)
    cpdef int getLastRefPos(self)
    cdef cython.str getCigarStringNP(self)
    cdef update_tagsNP(self)
    cpdef update_tags(self)

cpdef Layout_t MergeLayoutsToLayout(Layout_t, Layout_t)