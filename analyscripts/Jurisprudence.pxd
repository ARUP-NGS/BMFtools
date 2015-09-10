from libc.stdint cimport *

cimport pysam.calignmentfile
cimport pysam.calignedsegment
cimport cython

ctypedef pysam.calignmentfile.AlignmentFile AlignmentFile_t
ctypedef pysam.calignedsegment.AlignedSegment AlignedSegment_t
ctypedef cython.str cystr
ctypedef double double_t