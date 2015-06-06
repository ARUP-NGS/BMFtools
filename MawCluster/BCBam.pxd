cimport cython
cimport pysam.calignmentfile
cimport pysam.cfaidx
from cython cimport bint
from numpy cimport ndarray
from utilBMF.HTSUtils cimport cystr
ctypedef pysam.calignmentfile.AlignedSegment cAlignedSegment

cdef BarcodeTagCOBam_(pysam.calignmentfile.AlignmentFile bam,
                      pysam.calignmentfile.AlignmentFile outbam)
cpdef BarcodeTagCOBam(cystr bam, cystr outbam=?)

cdef dict GetCOTagDict_(cAlignedSegment read)

cpdef dict GetCOTagDict(cAlignedSegment read)

cdef double getAF(cAlignedSegment read)
cdef double getSF(cAlignedSegment read)
cdef cAlignedSegment TagAlignedSegment(cAlignedSegment read)