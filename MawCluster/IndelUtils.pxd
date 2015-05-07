cimport cython

cdef class IDVCFLine(object):

    """
    Contains data for writing to a VCF Line, where each ALT has its own line.

    minNumSS is the minumum number of start/stop combinations required to
    support a variant call.
    """
    
    cdef public cython.str TYPE, REF, ALT, ID, CHROM, FILTER, FormatStr, str
    cdef public cython.long POS, LEN, NumStartStops, NDPS, DPA
    cdef public cython.float reverseStrandFraction, QUAL, MDP
    cdef public cython.bint BothStrandSupport
    cdef public dict InfoFields, FormatFields
