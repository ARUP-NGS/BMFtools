cimport cython
cimport numpy as np
cimport pysam.calignmentfile
cimport pysam.cfaidx
cimport utilBMF.HTSUtils
from numpy cimport ndarray
from utilBMF.HTSUtils cimport cystr, chr2ph, TagTypeDict, BamTag, cReadsOverlap
from utilBMF.PysamUtils cimport CopyAlignedSegment
from utilBMF.cstring cimport cs_to_ia
from cython cimport bint
from cpython.array cimport array as py_array
from libc.stdlib cimport malloc, free, realloc
from libc.stdint cimport uint16_t
from utilBMF.Inliners cimport ph2chrInline, chrInline, CigarStrInline, opLenToStr
from utilBMF.PysamUtils cimport PysamToChrInline
from MawCluster.Probability cimport c_abs
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
    char MergeOverlappedLayouts(ArrayLayout_t AL1, ArrayLayout_t AL2)


cdef class Layout:
    cdef ArrayLayout_t Layout

    cdef public uchar mapq
    cdef public int tlen, pnext, flag, InitPos, firstMapped, reference_id, pos, rnext
    cdef public dict tagDict
    cdef public cystr Name
    cdef public bint is_reverse, isMerged, MergeSuccess

    cdef int getFirstAlignedRefPos(self)
    cdef int getLastAlignedRefPos(self)
    cdef int getFirstMappedReadPos(self)
    cdef int getAlignmentStart(self)

    cdef bint cPosIsMapped(self, int position)
    cpdef bint posIsMapped(self, int position)

    cdef char MergeLayouts_in_place(self, ArrayLayout_t pairedLayout)
    cdef ndarray[np.int16_t, ndim=1] getMergedPositions(self)
    cpdef ndarray[np.int16_t, ndim=1] getAgreement(self)
    cdef ndarray[np.int16_t, ndim=1] cGetAgreement(self)

    cpdef cystr getPVString(self)
    cdef cystr cGetPVString(self)
    cpdef cystr getFAString(self)
    cdef cystr cGetFAString(self)

    # Utilities for laying out.
    cpdef ndarray[np.int8_t, ndim=1] getOperations(self)
    cdef ndarray[np.int8_t, ndim=1] cGetOperations(self)

    # Utilities for creating output objects
    cdef cystr cGetCigarString(self)
    cpdef cystr getCigarString(self)

    cdef ndarray[int, ndim=1] cGetQual(self)
    cpdef ndarray[int, ndim=1] getQual(self)

    cdef cystr cGetQualString(self)
    cpdef cystr getQualString(self)

    cdef ndarray[char] cGetSeqArr(self)
    cpdef ndarray[char] getSeqArr(self)
    cpdef cystr getSeq(self)

    cdef list get_merge_tags(self)
    cdef get_tag_string(self)

    # Utilities for BMF metadata
    cdef ndarray[int, ndim=1] cGetGenomicDiscordantPositions(self)
    cdef ndarray[np.int16_t, ndim=1] cGetReadDiscordantPositions(self)
    cdef ndarray[np.int16_t, ndim=1] getMergeAgreements(self)

    # Output objects
    cpdef cystr __sam__(self, Layout_t pairedLayout)
    cdef cystr __csam__(self, Layout_t pairedLayout)

    # Updates
    cdef update_read_positions(self)
    cdef set_merge_tags_BT(self)

    # Slices
    cdef ndarray[np.int16_t, ndim=1] cGetAgreementSlice(self, size_t start=?)
    cdef ndarray[int, ndim=1] cGetQualSlice(self, size_t start=?)
    cdef cystr cGetCigarStringSlice(self, size_t start=?)
    cdef ndarray[char] cGetSeqArrSlice(self, size_t start=?)
    cdef ndarray[np.int8_t, ndim=1] cGetOperationsSlice(self, size_t start=?)
    cdef cystr cGetQualStringSlice(self, size_t start=?)


cdef public dict CigarDict, CigarStrDict
cdef int getLayoutLen(AlignedSegment_t read)  # Gets layout length
cpdef cystr FlattenCigarString(cystr cigar)
cdef cystr cFlattenCigarString(cystr cigar)
cpdef cystr ALPToStr(ArrayLayoutPos_t ALP)

cdef class ListBool:
    cdef list List
    cdef bint Bool


cdef object oagtag, oig0, oig1

cdef inline bint LayoutsOverlap(Layout_t Layout1, Layout_t Layout2):
    cdef int end1, end2, start1, start2
    end1 = Layout1.getLastAlignedRefPos()
    end2 = Layout2.getLastAlignedRefPos()
    start1 = Layout1.getFirstAlignedRefPos()
    start2 = Layout2.getFirstAlignedRefPos()
    if(start2 > end1 or start1 > end2):
        return False
    else:
        return True