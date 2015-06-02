cimport cython
cimport numpy as np
cimport pysam.calignmentfile
from utilBMF.HTSUtils cimport pPileupRead
ctypedef pPileupColumn pPileupColumn_t
ctypedef pPileupRead pPileupRead_t


cdef class AlleleAggregateInfo:

    """
    Class which holds summary information for a given alt allele
    at a specific position. Meant to be used as part of the PCInfo
    class as values in the VariantDict.
    All alt alleles in this set of reads should be identical.
    recList must be a list of PRInfo objects.

    """
    cdef public cython.long DOC
    cdef public cython.long DOCTotal
    cdef public cython.long NUMALT
    cdef public cython.long pos
    cdef public cython.long len
    cdef public cython.long TND
    cdef public cython.long TotalReads
    cdef public cython.long MergedReads
    cdef public cython.long ReverseMergedReads
    cdef public cython.long ForwardMergedReads
    cdef public cython.long ReverseTotalReads
    cdef public cython.long ForwardTotalReads
    cdef public cython.long SumBQScore
    cdef public cython.long SumMQScore
    cdef public cython.long FSR
    cdef public cython.long NumberDuplexReads, maxND
    cdef public cython.float MNF, maxNF, NFSD, AveFamSize, AveMQ, AveBQ, minMQ, minBQ
    cdef public cython.float minFA, MFractionAgreed, reverseStrandFraction, MFA
    cdef public cython.float minFrac, TotalAlleleFrequency, MergedAlleleFrequency
    cdef public cython.float minPVFrac, AABPSD, AAMBP, MBP, BPSD, PFSD, MPF
    cdef public list recList
    cdef public dict TotalAlleleDict, StrandCountsDict, StrandCountsTotalDict
    cdef public dict strandedTransitionDict
    cdef public cython.str ALT, consensus, transition, contig
    cdef public cython.bint BothStrandSupport


cdef class pPileupColumn:
    """
    Python container for the PileupColumn proxy in pysam.
    """
    cdef public cython.long nsegments
    cdef public cython.long reference_id
    cdef public cython.long reference_pos
    cdef public list pileups


cdef class PRInfo:

    """
    Created from a pysam.PileupRead object or its python equivalent,
    a pPileupRead.
    Holds family size, SV tags, base quality, mapping quality, and base.
    If any of the "finicky" fields aren't filled (e.g., if BAMs are not
    produced using BMFTools), they are set to None. Check to see if an attribute
    is None first if you'd like to save yourself some headaches.
    """
    cdef public cython.bint Pass, is_reverse, is_proper_pair
    cdef public cython.long FM, BQ, MQ, query_position, FA, PV, ND
    cdef public pysam.calignmentfile.AlignedSegment read
    cdef public cython.str SVTagString
    cdef public cython.str BaseCall
    cdef public cython.str ssString, query_name
    cdef public cython.float FractionAgreed, PVFrac, NF
    cdef public np.ndarray PV_Array
    cpdef object opt(self, cython.str arg)


cdef class PCInfo:

    """
    Takes a pysam.PileupColumn covering one base in the reference
    and makes a new class which has "reference" (inferred from
    consensus) and a list of PRData (one for each read).
    The option "duplex required" determines whether or not variants must
    be supported by both reads in a pair for a proper call. Increases
    stringency, might lose some sensitivity for a given sequencing depth.

    exclusionSVTags should be a string of comma-separated SV tags.
    The presence of one of these tags in a read causes it to be thrown out
    of the pileup.
    """
    cdef public cython.str experiment, excludedSVTagStr, consensus, TotalFracStr
    cdef public cython.str MergedFracStr, MergedCountStr, TotalCountStr, str
    cdef public cython.str MergedStrandednessStr, TotalStrandednessStr, AlleleFreqStr
    cdef public cython.str contig
    cdef public cython.long minMQ, minBQ, pos, FailedQCReads, FailedFMReads
    cdef public cython.long FailedBQReads, FailedAFReads, FailedMQReads
    cdef public cython.long FailedSVReads, MergedReads, TotalReads
    cdef public cython.float minAF, reverseStrandFraction, AAMBP, AABPSD
    cdef public pPileupColumn_t PCol
    cdef public dict FailedSVReadDict, VariantDict, TotalAlleleFreqDict
    cdef public dict MergedAlleleFreqDict, MergedAlleleDict, TotalAlleleDict
    cdef public dict MergedStrandednessRatioDict, TotalStrandednessRatioDict
    cdef public list Records, AltAlleleData, DiscNames
    cdef public cython.bint BothStrandAlignment
    cdef public cython.long maxND, FailedNDReads
