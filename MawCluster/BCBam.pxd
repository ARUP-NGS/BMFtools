cimport cython
cimport pysam.calignmentfile
cimport pysam.cfaidx
from cython cimport bint
from numpy cimport ndarray
ctypedef pysam.calignmentfile.AlignedSegment cAlignedSegment

cdef BarcodeTagCOBam_(pysam.calignmentfile.AlignmentFile bam,
                      pysam.calignmentfile.AlignmentFile outbam,
                      bint addRG=?)
cpdef BarcodeTagCOBam(cython.str bam, cython.str realigner=?)

cdef dict GetCOTagDict_(cAlignedSegment read)

cpdef dict GetCOTagDict(cAlignedSegment read)

cdef double getAF(cAlignedSegment read)
cdef double getSF(cAlignedSegment read)