cimport cython
cimport pysam.calignmentfile
cimport pysam.cfaidx
from cython cimport bint
from numpy cimport ndarray
from utilBMF.HTSUtils cimport PysamToChrDict
from utilBMF.Inliners cimport Num2Nuc
from utilBMF.PysamUtils cimport PysamToChrInline
ctypedef pysam.calignmentfile.AlignedSegment AlignedSegment_t
ctypedef pysam.calignmentfile.AlignmentFile AlignmentFile
ctypedef pysam.calignmentfile.AlignmentFile AlignmentFile_t
ctypedef cython.str cystr

cdef cystr cBarcodeTagCOBam(pysam.calignmentfile.AlignmentFile bam,
                      pysam.calignmentfile.AlignmentFile outbam)
cpdef AlignedSegment_t pTagAlignedSegment(AlignedSegment_t read)
cpdef cystr pBarcodeTagCOBam(cystr bam, cystr outbam=?)

cdef dict cGetCOTagDict(AlignedSegment_t read)

cpdef dict pGetCOTagDict(AlignedSegment_t read)

cdef double getAF(AlignedSegment_t read)
cdef double getSF(AlignedSegment_t read)
# cpdef AlignedSegment_t TagAlignedSegment(AlignedSegment_t read)


cdef inline cystr RPString(AlignedSegment_t read):
    return (PysamToChrInline(read.reference_id) + ":%s," % read.pos +
            PysamToChrInline(read.rnext) +
            ":%s" % read.mpos)

cdef class BamPipe:
    """
    Creates a callable function which acts on a BAM stream.

    :param function - callable function which returns an input BAM object.
    :param bin_input - boolean - true if input is BAM
    false for TAM/SAM
    :param bin_output - boolean - true to output in BAM format.
    :param uncompressed_output - boolean - true to output uncompressed
    BAM records.
    """
    cdef public object function
    cdef public bint bin_input, bin_output, uncompressed_output
    cdef public AlignmentFile_t inHandle, outHandle
    cpdef process(self)