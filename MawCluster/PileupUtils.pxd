cimport cython
from cython cimport bint
cimport numpy as np
cimport pysam.calignmentfile
cimport utilBMF.HTSUtils
from numpy cimport ndarray
from utilBMF.HTSUtils cimport cystr, nucList, PysamToChrDict
from utilBMF.HTSUtils cimport PileupReadPair
from utilBMF.HTSUtils cimport pPileupRead
from utilBMF.PysamUtils cimport PysamToChrInline
from cpython cimport array as c_array
ctypedef np.int32_t np_int32_t
ctypedef c_array.array py_array
ctypedef AlleleAggregateInfo AlleleAggregateInfo_t
ctypedef PileupReadPair PileupReadPair_t
ctypedef pPileupColumn pPileupColumn_t
ctypedef pPileupRead pPileupRead_t
ctypedef PRInfo PRInfo_t
ctypedef pysam.calignmentfile.PileupRead PileupRead_t
ctypedef pysam.calignmentfile.AlignedSegment AlignedSegment_t
ctypedef pysam.calignmentfile.AlignmentFile AlignmentFile_t


cdef class AlleleAggregateInfo:

    """
    Class which holds summary information for a given alt allele
    at a specific position. Meant to be used as part of the PCInfo
    class as values in the VariantDict.
    All alt alleles in this set of reads should be identical.
    recList must be a list of PRInfo objects.

    """
    cdef public long DOC
    cdef public long DOCTotal
    cdef public long NUMALT
    cdef public long pos
    cdef public long len
    cdef public long TND
    cdef public long TotalReads
    cdef public long MergedReads
    cdef public long ReverseMergedReads
    cdef public long ForwardMergedReads
    cdef public long ReverseTotalReads
    cdef public long ForwardTotalReads
    cdef public long SumBQScore
    cdef public long SumMQScore
    cdef public long FSR
    cdef public long NumberDuplexReads, maxND
    cdef public cython.float MNF, maxNF, NFSD, AveFamSize, AveMQ, AveBQ, minMQ, minBQ
    cdef public cython.float minFA, MFractionAgreed, reverseStrandFraction, MFA
    cdef public cython.float minFrac, TotalAlleleFrequency, MergedAlleleFrequency
    cdef public cython.float minPVFrac, AABPSD, AAMBP, MBP, BPSD, PFSD, MPF
    cdef public list recList
    cdef public dict TotalAlleleDict, StrandCountsDict, StrandCountsTotalDict
    cdef public dict strandedTransitionDict
    cdef public cystr ALT, consensus, transition, contig
    cdef public bint BothStrandSupport


cdef class pPileupColumn:
    """
    Python container for the PileupColumn proxy in pysam.
    """
    cdef public long nsegments
    cdef public long reference_id
    cdef public long reference_pos
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
    cdef public bint Pass, is_reverse, is_proper_pair
    cdef public bint SVPass
    cdef public long FM, MQ, query_position, FA, ND
    cdef public np_int32_t BQ
    cdef public AlignedSegment_t read
    cdef public cystr BaseCall
    cdef public cystr ssString, query_name
    cdef public cython.float FractionAgreed, PVFrac, NF, AF
    cpdef object opt(self, cystr arg)


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
    cdef public cystr experiment, excludedSVTagStr, consensus, TotalFracStr
    cdef public cystr MergedFracStr, MergedCountStr, TotalCountStr, str
    cdef public cystr MergedStrandednessStr, TotalStrandednessStr, AlleleFreqStr
    cdef public cystr contig
    cdef public long minMQ, minBQ, pos, FailedQCReads, FailedFAReads
    cdef public long FailedBQReads, FailedAFReads, FailedMQReads, FailedNDReads
    cdef public long FailedAMPReads
    cdef public long FailedSVReads, MergedReads, TotalReads
    cdef public cython.float minAF, reverseStrandFraction, AAMBP, AABPSD
    cdef public pPileupColumn_t PCol
    cdef public dict FailedSVReadDict, VariantDict, TotalAlleleFreqDict
    cdef public dict MergedAlleleFreqDict, MergedAlleleDict, TotalAlleleDict
    cdef public dict MergedStrandednessRatioDict, TotalStrandednessRatioDict
    cdef public list Records, AltAlleleData, DiscNames
    cdef public bint BothStrandAlignment
    cdef public long maxND

cdef inline bint TestSVFlag(cystr tagflag):
    if(tagflag == "LI"):
        return False
    elif(tagflag == "MI"):
        return False
    elif(tagflag == "MDC"):
        return False
    elif(tagflag == "MSS"):
        return False
    elif(tagflag == "ORU"):
        return True
    elif(tagflag == "ORB"):
        return True
    elif(tagflag == "ORS"):
        return True
    elif(tagflag == "DRP"):
        return True
    elif(tagflag == "DSI"):
        return True
    elif(tagflag == "DSD"):
        return True
    elif(tagflag == "NF"):
        return True
    else:
        return True

cdef inline bint TestSVTags(cystr inTag):
    cdef cystr tagflag
    for tagflag in inTag.split(","):
        if(TestSVFlag(tagflag) is False):
            return False
    return True


cdef inline cystr AS_to_key(AlignedSegment_t read):
    if(read.is_read1):
        return read.query_name + "1"
    else:
        return read.query_name + "2"