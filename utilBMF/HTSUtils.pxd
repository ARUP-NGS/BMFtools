cimport cython
cimport pysam.calignmentfile
cimport numpy as np
ctypedef PileupReadPair PileupReadPair_t
ctypedef np.longdouble_t dtype128_t
ctypedef pPileupRead pPileupRead_t
ctypedef ReadPair ReadPair_t

cdef class pPileupRead:
    """
    Python container for the PileupRead proxy in pysam
    """
    cdef public cython.str BaseCall
    cdef public cython.bint is_del
    cdef public cython.long level
    cdef public cython.long indel
    cdef public cython.long query_position
    cdef public cython.str name, str
    cdef public pysam.calignmentfile.AlignedSegment alignment

cdef class PileupReadPair:

    """
    Holds both bam record objects in a pair of pileup reads.
    Currently, one read unmapped and one read soft-clipped are
    both marked as soft-clipped reads.
    Accepts a list of length two as input.
    """

    cdef public pPileupRead_t read1
    cdef public pPileupRead_t read2
    cdef public ReadPair_t RP
    cdef public cython.bint discordant
    cdef public cython.str discordanceString
    cdef public cython.str name


cdef class ReadPair:

    """
    Holds both bam record objects in a pair.
    Currently, one read unmapped and one read soft-clipped are
    both marked as soft-clipped reads.
    """
    cdef public pysam.calignmentfile.AlignedSegment read1
    cdef public pysam.calignmentfile.AlignedSegment read2
    cdef public list SVTags
    cdef public cython.bint read1_is_unmapped
    cdef public cython.bint read1_soft_clipped
    cdef public cython.bint read2_is_unmapped
    cdef public cython.bint read2_soft_clipped
    cdef public cython.bint SameContig
    cdef public cython.str read1_contig
    cdef public cython.str read2_contig
    cdef public cython.str ContigString
    cdef public cython.long insert_size
    cdef public cython.bint read1_in_bed
    cdef public cython.bint read2_in_bed


cdef class AbstractIndelContainer:
    """
    Base class for insertion and deletion container objects.

    Type can be -1, 0, or 1. (Deletion, deletion and insertion, and just insertion)
    Start and end refer to different things for insertions and deletions.
    For a deletion, start is the first missing base and end is the last missing
    reference base position.
    seq should be None for a deletion
    """
    cdef public cython.str contig
    cdef public cython.long start
    cdef public cython.long end
    cdef public cython.long type
    cdef public cython.str seq
    cdef public cython.float shen
    cdef public list readnames
    cdef public cython.str uniqStr

cdef class Deletion(AbstractIndelContainer):
    """
    See HTSUtils.pyx for doc string.
    """


cdef class Insertion(AbstractIndelContainer):
    """
    See HTSUtils.pyx for doc string.
    """
