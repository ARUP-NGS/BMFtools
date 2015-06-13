cimport cython
cimport pysam.calignmentfile
cimport pysam.cfaidx
from cython cimport bint
from numpy cimport ndarray
from utilBMF.HTSUtils cimport cystr
from utilBMF.HTSUtils cimport PysamToChrDict
from MawCluster.BCFastq cimport letterNumDict
ctypedef pysam.calignmentfile.AlignedSegment cAlignedSegment
ctypedef pysam.calignmentfile.AlignmentFile AlignmentFile

cdef cystr cBarcodeTagCOBam(pysam.calignmentfile.AlignmentFile bam,
                      pysam.calignmentfile.AlignmentFile outbam)
cpdef cystr pBarcodeTagCOBam(cystr bam, cystr outbam=?)

cdef dict cGetCOTagDict(cAlignedSegment read)

cpdef dict pGetCOTagDict(cAlignedSegment read)

cdef double getAF(cAlignedSegment read)
cdef double getSF(cAlignedSegment read)
# cpdef cAlignedSegment TagAlignedSegment(cAlignedSegment read)