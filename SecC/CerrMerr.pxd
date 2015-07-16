cimport cython
cimport numpy as np
from libc.stdint cimport int8_t, int32_t
from utilBMF.Inliners cimport CONTEXT_TO_ARRAY_POS

cimport pysam.calignedsegment
ctypedef pysam.calignedsegment.AlignedSegment AlignedSegment_t
ctypedef cython.str cystr
from numpy cimport ndarray
