cimport cython
cimport numpy as np
ctypedef np.longdouble_t dtype128_t


cdef class SNVCFLine:

    """
    Contains data for writing to a VCF Line, where each ALT has its own line.

    If reverseStrandFraction is negative, it means that that information was
    not provided at initialization.
    minNumSS is the minumum number of start/stop combinations required to
    support a variant call.
    """
    cdef public cython.str REF, CHROM, CONS, ALT, ID, FILTER, InfoStr, FormatStr, str
    cdef public cython.str FormatKey, FormatValue
    cdef public cython.long NumStartStops, POS
    cdef public cython.bint BothStrandSupport, AABothStrandAlignment
    cdef public dtype128_t reverseStrandFraction, QUAL
    cdef public dict InfoFields, FormatFields



cdef class VCFPos:
    """
    Holds the multiple SNVCFLine Objects and other information
    for writing VCF lines for each alt at a given position.
    """
    cdef public cython.long pos, minMQ, FailedAFReads, minDuplexPairs
    cdef public cython.str consensus, TotalFracStr, MergedFracStr, TotalCountStr
    cdef public cython.str MergedCountStr, REF, EST, str
    cdef public cython.float reverseStrandFraction, minAF
    cdef public cython.bint AABothStrandAlignment, requireDuplex, keepConsensus
    cdef public cython.bint DuplexRequired
    cdef public list VCFLines


