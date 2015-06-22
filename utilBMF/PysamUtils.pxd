cimport pysam.calignmentfile
from cpython cimport array
ctypedef array.array py_array

ctypedef pysam.calignmentfile.AlignedSegment AlignedSegment_t
ctypedef pysam.calignmentfile.AlignmentFile AlignmentFile_t

cdef AlignedSegment_t CopyAlignedSegment(AlignedSegment_t template)

cdef class PairwiseAlignmentFile:
    cdef public AlignmentFile_t handle