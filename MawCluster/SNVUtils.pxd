cimport cython
cimport numpy as np
cimport pysam.cfaidx
cimport pysam.calignmentfile
ctypedef np.longdouble_t dtype128_t
ctypedef cython.str cystr
ctypedef SNVCFLine SNVCFLine_t
from MawCluster.PileupUtils cimport AlleleAggregateInfo, PCInfo
ctypedef AlleleAggregateInfo AlleleAggregateInfo_t
ctypedef PCInfo PCInfo_t
ctypedef pysam.calignmentfile.PileupRead PileupRead


cdef class SNVCFLine:

    """
    Contains data for writing to a VCF Line, where each ALT has its own line.

    If reverseStrandFraction is negative, it means that that information was
    not provided at initialization.
    minNumSS is the minumum number of start/stop combinations required to
    support a variant call.
    """
    cdef public cystr REF, CHROM, CONS, ALT, ID, FILTER, InfoStr, FormatStr
    cdef public cystr str
    cdef public cystr FormatKey, FormatValue
    cdef public cystr flankingBefore, flankingAfter
    cdef public long NumStartStops, POS
    cdef public cython.bint BothStrandSupport, AABothStrandAlignment
    cdef public dtype128_t reverseStrandFraction, QUAL
    cdef public dict InfoFields, FormatFields


cdef class VCFPos:
    """
    Holds the multiple SNVCFLine Objects and other information
    for writing VCF lines for each alt at a given position.
    """
    cdef public long pos, minMQ, FailedAFReads, minDuplexPairs
    cdef public cystr consensus, TotalFracStr, MergedFracStr, TotalCountStr
    cdef public cystr MergedCountStr, REF, EST, str
    cdef public cython.float reverseStrandFraction, minAF
    cdef public cython.bint AABothStrandAlignment, requireDuplex, keepConsensus
    cdef public cython.bint DuplexRequired
    cdef public list VCFLines
