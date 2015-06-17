cimport cython
cimport numpy as np
cimport pysam.calignmentfile
cimport pysam.cfaidx
cimport utilBMF.HTSUtils
from numpy cimport ndarray
from utilBMF.HTSUtils cimport cystr, chr2ph, ph2chrDict, PysamToChrDict, TagTypeDict, BamTag
from cython cimport bint
from libc.stdlib cimport malloc, free, realloc
from libc.stdint cimport uint16_t
ctypedef pysam.calignmentfile.AlignedSegment AlignedSegment_t
ctypedef unsigned char uchar
ctypedef BamTag BamTag_t
ctypedef Layout Layout_t
ctypedef ArrayLayoutPos ArrayLayoutPos_t
ctypedef ArrayLayout ArrayLayout_t
ctypedef MergeRet MergeRet_t


cdef extern from "MPA.h":
    struct ArrayLayoutPos:
        int pos
        uint16_t readPos
        int quality
        uint16_t agreement
        char operation
        char base
        char mergeAgreed
    struct ArrayLayout:
        ArrayLayoutPos_t * layouts
        size_t length
    struct MergeRet:
        ArrayLayout_t Layout
        char Success
    ArrayLayoutPos_t cMergeLayoutPositions(ArrayLayoutPos_t L1, ArrayLayoutPos_t L2)
    int getFirstAlignedRefPos(ArrayLayout_t layout)
    ArrayLayout_t MergeLayouts(ArrayLayout_t AL1, ArrayLayout_t AL2)
    MergeRet_t MergeWithPassFail(ArrayLayout_t AL1, ArrayLayout_t AL2)


cdef class Layout:
    cdef ArrayLayout_t Layout

    cdef public uchar mapq
    cdef public int tlen, pnext, flag, InitPos, firstMapped, reference_id, pos, rnext
    cdef public dict tagDict
    cdef public cystr Name
    cdef public bint is_reverse, isMerged, MergeSuccess

    cdef int getFirstAlignedRefPos(self)
    cdef int getFirstMappedReadPos(self)
    cdef int getAlignmentStart(self)

    cdef bint cPosIsMapped(self, int position)
    cpdef bint posIsMapped(self, int position)

    cdef MergeLayouts_in_place(self, ArrayLayout_t pairedLayout)
    cpdef MergeLayout(self, Layout_t pairedLayout)
    cdef update_tags_(self)
    cdef ndarray[np.int16_t, ndim=1] getMergedPositions(self)
    cpdef ndarray[int, ndim=1] getAgreement(self)
    cdef ndarray[uint16_t, ndim=1] cGetAgreement(self)

    # Utilities for laying out.
    cpdef ndarray[np.int8_t, ndim=1] getOperations(self)
    cdef ndarray[np.int8_t, ndim=1] cGetOperations(self)

    # Utilities for creating output objects
    cdef cystr cGetCigarString(self)
    cpdef cystr getCigarString(self)
    cdef cystr cGetCigarString2(self)
    cpdef cystr getCigarString2(self)

    cdef ndarray[np.int16_t, ndim=1] cGetQual(self)
    cpdef ndarray[int, ndim=1] getQual(self)

    cdef cystr cGetQualString(self)
    cpdef cystr getQualString(self)

    cdef ndarray[char] cGetSeqArr(self)
    cpdef ndarray[char] getSeqArr(self)
    cpdef cystr getSeq(self)

    # Utilities for BMF metadata
    cdef ndarray[int, ndim=1] cGetGenomicDiscordantPositions(self)
    cdef ndarray[int, ndim=1] cGetReadDiscordantPositions(self)
    cdef ndarray[np.int16_t, ndim=1] getMergeAgreements(self)

    # Output objects
    cdef AlignedSegment_t __read__(self)

    # Updates
    cdef update_read_positions(self)


cdef public dict chrDict, CigarDict, CigarStrDict
cdef int getLayoutLen(AlignedSegment_t read)  # Gets layout length
cpdef cystr ALPToStr(ArrayLayoutPos_t ALP)

cdef class ListBool:
    cdef list List
    cdef bint Bool

cdef object oagtag