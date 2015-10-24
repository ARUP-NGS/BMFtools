cimport cython
cimport pysam.calignmentfile
cimport pysam.calignedsegment
cimport pysam.cfaidx
cimport numpy as np
from numpy cimport ndarray
from cython cimport bint
from utilBMF.cstring cimport cs_to_ph, cs_to_ia
from cpython cimport array as c_array
from pysam.cfaidx cimport PersistentFastqProxy
from libc.stdint cimport int8_t, uint32_t, uint64_t
from libc.stdlib cimport *
from libc.stdio cimport *
ctypedef c_array.array py_array
ctypedef cython.str cystr
ctypedef np.longdouble_t dtype128_t
ctypedef pPileupRead pPileupRead_t
ctypedef pysam.calignedsegment.AlignedSegment AlignedSegment_t
ctypedef pFastqProxy pFastqProxy_t
ctypedef pysam.cfaidx.FastqProxy FastqProxy_t

cimport pysam.ctabixproxies
ctypedef pysam.calignedsegment.PileupRead cPileupRead
ctypedef Insertion Insertion_t
ctypedef Deletion Deletion_t
ctypedef AbstractIndelContainer AbstractIndelContainer_t
ctypedef IndelQuiver IndelQuiver_t
ctypedef IDVCFLine IDVCFLine_t
ctypedef np.int32_t np_int32_t

makeFastqProxy = pysam.cfaidx.makeFastqProxy

cdef class pPileupRead:
    """
    Python container for the PileupRead proxy in pysam
    """
    cdef public cystr BaseCall
    cdef public bint is_del
    cdef public long level
    cdef public long indel
    cdef public long query_position
    cdef public cystr name, str
    cdef public AlignedSegment_t alignment
    cdef public np_int32_t BQ, FA, MBQ, FM, MQ
    cpdef object opt(self, cystr arg)


cdef class AbstractIndelContainer:
    """
    Base class for insertion and deletion container objects.

    Type can be -1, 0, or 1. (Deletion, deletion and insertion, and
    just insertion)
    Start and end refer to different things for insertions and deletions.
    For a deletion, start is the first missing base and end is the last
    missing reference base position.

    seq should be None for a deletion
    """
    cdef public cystr contig, seq, uniqStr
    cdef public long type, shenwindow, end, start
    cdef public cython.float shen
    cdef public list readnames, StartStops
    cdef public pysam.cfaidx.FastaFile handle

cdef class Deletion(AbstractIndelContainer):
    """
    See HTSUtils.pyx for doc string.
    """


cdef class Insertion(AbstractIndelContainer):
    """
    See HTSUtils.pyx for doc string.
    """

cdef class IndelQuiver(object):
    cdef public dict data, counts, readnames
    cdef public long window, minMQ, minFM, minPairs, minNumSS
    cdef public pysam.cfaidx.FastaFile fastaRef
    cdef public pysam.calignmentfile.AlignmentFile bam
    cdef public cython.float minShen


cdef class IDVCFLine(object):

    """
    Contains data for writing to a VCF Line, where each ALT has its own line.

    minNumSS is the minumum number of start/stop combinations required to
    support a variant call.
    """

    cdef public cystr TYPE, REF, ALT, ID, CHROM, FILTER, FormatStr, str
    cdef public long POS, LEN, NumStartStops, NDPS, DPA
    cdef public cython.float reverseStrandFraction, QUAL, MDP
    cdef public bint BothStrandSupport
    cdef public dict InfoFields, FormatFields


cdef class pFastqFile(object):
    cdef public pysam.cfaidx.FastqFile handle
    cpdef close(self)

cdef class BamTag(object):
    """
    Contains a tag, a value, and a type, all of which are string objects.
    """
    cdef readonly cystr tag
    cdef readonly cystr tagtype
    cdef public object value

cpdef public cystr RevCmp(cystr seq)

cpdef public list permuteNucleotides(long maxn, object nci=?, int kmerLen=?)

cpdef bint ReadsOverlap(
        AlignedSegment_t read1,
        AlignedSegment_t read2)


cdef bint cReadsOverlap(
        AlignedSegment_t read1,
        AlignedSegment_t read2)

cdef double cyOptStdDev_(ndarray[np.float64_t, ndim=1] a)
cdef public cystr cGetBS(pFastqProxy_t)


cdef public dict PysamToChrDict, ph2chrDict, TypeConversionDict
cdef public dict chr2ph, chr2phStr, int2Str, TagTypeDict
cdef public list nucList

cdef extern from "../lib/cstr_utils.h" nogil:
    char *revcmp(char *, uint64_t);

cdef class pFastqProxy:
    cdef public cystr comment, name, quality, sequence
    cdef cystr tostring(self)
    cpdef dict get_tag_dict(self)
    cpdef int get_int_tag(self, cystr tag)
    cpdef py_array getQualArray(self)
    cpdef cystr getSlice(self, int start=?, int end=?,
                         cystr addComment=?)
