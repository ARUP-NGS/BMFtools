cimport cython
cimport pysam.calignmentfile
ctypedef pysam.calignmentfile.AlignedSegment cAlignedSegment

cdef BarcodeTagCOBam_(pysam.calignmentfile.AlignmentFile bam,
                      pysam.calignmentfile.AlignmentFile outbam,
                      cython.bint addRG=?)
cpdef BarcodeTagCOBam(cython.str bam, cython.str realigner=?)

cdef dict GetCOTagDict_(cAlignedSegment read)

cpdef dict GetCOTagDict(cAlignedSegment read)