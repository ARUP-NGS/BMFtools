cimport cython
cimport numpy as np
from libc.stdint cimport int8_t, int32_t, int64_t, int16_t
from utilBMF.Inliners cimport CONTEXT_TO_ARRAY_POS

cimport pysam.calignedsegment
cimport pysam.calignmentfile
ctypedef pysam.calignedsegment.AlignedSegment AlignedSegment_t
ctypedef pysam.calignmentfile.AlignmentFile AlignmentFile_t
ctypedef cython.str cystr
ctypedef double double_t
from numpy cimport ndarray
