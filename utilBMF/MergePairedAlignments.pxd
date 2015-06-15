cimport cython
cimport numpy as np
cimport pysam.calignmentfile
cimport pysam.cfaidx
cimport utilBMF.HTSUtils
from numpy cimport ndarray
from utilBMF.HTSUtils cimport cystr, chr2ph, ph2chrDict, PysamToChrDict, TagTypeDict, cReadsOverlap, BamTag
from cython cimport bint
from libc.stdlib cimport malloc, free
ctypedef Layout Layout_t
ctypedef LayoutPos LayoutPos_t
ctypedef ArrayLayout ArrayLayout_t
ctypedef pysam.calignmentfile.AlignedSegment AlignedSegment_t
ctypedef unsigned char uchar
ctypedef BamTag BamTag_t

cdef struct ArrayLayoutPos:
    int pos
    int readPos
    int quality
    int agreement
    char operation
    char base
    char mergeAgreed
ctypedef ArrayLayoutPos ArrayLayoutPos_t

cdef class LayoutPos:
    cdef public cython.int pos, readPos, quality, agreement
    cdef public char operation, base, mergeAgreed
    cdef bint merged
    # cdef public cystr operation, base
    cpdef bint ismapped(self)
    cdef bint getMergeAgreed(self)
    cdef bint getMergeSet(self)

cdef class ArrayLayout:
    cdef public AlignedSegment_t read
    cdef ArrayLayoutPos_t *layoutArray
    cdef public size_t length

    cdef public uchar mapq
    cdef public int tlen, pnext, flag, InitPos, firstMapped

    cdef public dict tagDict

    cdef bint cPosIsMapped(self, int position)
    cpdef bint posIsMapped(self, int position)
    cdef int getFirstMappedRefPos(self)
    cdef ndarray[int, ndim=1] cGetQual(self)
    cpdef ndarray[int, ndim=1] getQual(self)
    cdef cystr cGetQualString(self)
    cpdef cystr getQualString(self)
    cdef ndarray[char] getSeqArr(self)
    cpdef cystr getSeq(self)

    
'''
cdef class ArrayLayoutPos:
    cdef public int[6] values
    # pos, readPos, quality, agreement, operation, base
    # cdef public cystr operation, base
    cpdef bint ismapped(self)

cdef class ArrayLayout:
    cdef public pysam.calignmentfile.AlignedSegment read
    cdef ndarray data
    cdef public bint isMerged, is_reverse
'''


cdef class Layout:
    cdef public pysam.calignmentfile.AlignedSegment read
    cdef public list positions
    cdef public dict tagDict
    cdef public cython.int firstMapped, InitPos, flag, pnext, tlen, mapq
    cdef public cystr Name, contig, rnext
    cdef public bint isMerged, is_reverse

    cpdef cython.int getAlignmentStart(self)
    cpdef cystr getCigarString(self)
    cpdef cystr getSeq(self)
    cdef ndarray[char] getSeqArr(self)
    cdef int cGetRefPosForFirstPos(self)
    cpdef int pGetRefPosForFirstPos(self)
    cpdef ndarray[int, ndim=1] getAgreement(self)
    cdef ndarray[int, ndim=1] cGetAgreement(self)
    cdef ndarray[int, ndim=1] cGetQual(self)
    cpdef ndarray[int, ndim=1] getQual(self)
    cdef cystr cGetQualString(self)
    cpdef cystr getQualString(self)
    cdef int cGetLastRefPos(self)
    cpdef int getLastRefPos(self)
    cdef cystr cGetCigarString(self)
    cdef update_tags_(self)
    cpdef update_tags(self)
    cdef ndarray[char, ndim=1] cGetMergedPositions(self)
    cdef ndarray[char, ndim=1] cGetMergeAgreements(self)
    cpdef ndarray[char, ndim=1] getMergedPositions(self)
    cpdef ndarray[char, ndim=1] getMergeAgreements(self)
    cdef ndarray[int, ndim=1] cGetGenomicDiscordantPositions(self)
    cdef ndarray[int, ndim=1] cGetReadDiscordantPositions(self)

cpdef Layout_t MergeLayoutsToLayout(Layout_t, Layout_t)
cdef list CigarOpToLayoutPosList(int offset, int cigarOp, int cigarLen,
                                 pysam.calignmentfile.AlignedSegment rec)

cdef public dict chrDict, CigarDict
cdef LayoutPos_t cMergePositions(LayoutPos_t pos1, LayoutPos_t pos2)
cdef int getLayoutLen(AlignedSegment_t read)

cdef class ListBool:
    cdef list List
    cdef bint Bool
