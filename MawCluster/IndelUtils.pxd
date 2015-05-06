cimport cython
cimport pysam.cfaidx
cimport pysam.calignmentfile

cdef class AbstractIndelContainer:
    """
    Base class for insertion and deletion container objects.

    Type can be -1, 0, or 1. (Deletion, deletion and insertion, and just insertion)
    Start and end refer to different things for insertions and deletions.
    For a deletion, start is the first missing base and end is the last missing
    reference base position.
    seq should be None for a deletion
    """
    cdef public cython.str contig, seq, uniqStr
    cdef public cython.long type, shenwindow, end, start
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
    cdef public cython.long window, minMQ, minFM, minPairs, minNumSS
    cdef public pysam.cfaidx.FastaFile fastaRef
    cdef public pysam.calignmentfile.AlignmentFile bam
    cdef public cython.float minShen