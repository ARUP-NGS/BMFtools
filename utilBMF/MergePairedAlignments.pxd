cimport cython
cimport numpy as np
cimport pysam.calignmentfile
cimport pysam.cfaidx
cimport utilBMF.HTSUtils
from numpy cimport ndarray
from utilBMF.HTSUtils cimport cystr, chr2ph, ph2chrDict, PysamToChrDict, TagTypeDict, BamTag
from cython cimport bint
from libc.stdlib cimport malloc, free, realloc
ctypedef pysam.calignmentfile.AlignedSegment AlignedSegment_t
ctypedef unsigned char uchar
ctypedef BamTag BamTag_t
ctypedef Layout Layout_t
ctypedef ArrayLayoutPos ArrayLayoutPos_t
ctypedef ArrayLayout ArrayLayout_t


cdef extern from "MPA.h":
    struct ArrayLayoutPos:
        int pos
        int readPos
        int quality
        int agreement
        char operation
        char base
        char mergeAgreed
    struct ArrayLayout:
        ArrayLayoutPos_t * layouts
        size_t length
    ArrayLayoutPos_t cMergeLayoutPositions(ArrayLayoutPos_t L1, ArrayLayoutPos_t L2)
    int getFirstAlignedRefPos(ArrayLayout_t layout)
    ArrayLayout_t MergeLayouts(ArrayLayout_t AL1, ArrayLayout_t AL2)


cdef class Layout:
    cdef ArrayLayout_t Layout

    cdef public uchar mapq
    cdef public int tlen, pnext, flag, InitPos, firstMapped
    cdef public dict tagDict

    cdef int getFirstAlignedRefPos(self)

    cdef bint cPosIsMapped(self, int position)
    cpdef bint posIsMapped(self, int position)
    cdef int getFirstMappedRefPos(self)
    cdef ndarray[np.int16_t, ndim=1] cGetQual(self)
    cpdef ndarray[int, ndim=1] getQual(self)
    cdef cystr cGetQualString(self)
    cpdef cystr getQualString(self)
    cdef ndarray[char] cGetSeqArr(self)
    cpdef ndarray[char] getSeqArr(self)
    cpdef cystr getSeq(self)
    cdef int getFirstMappedReadPos(self)
    cdef int getFirstMappedRefPos(self)
    cdef MergeLayouts_in_place(self, ArrayLayout_t pairedLayout)

cdef public dict chrDict, CigarDict
cdef int getLayoutLen(AlignedSegment_t read)  # Gets layout length
cpdef cystr ALPToStr(ArrayLayoutPos_t ALP)

cdef class ListBool:
    cdef list List
    cdef bint Bool
