cimport utilBMF.HTSUtils
cimport cython
cimport pysam.calignmentfile

from utilBMF.HTSUtils cimport cGetBS as GetFqBS
ctypedef pysam.calignmentfile.AlignedSegment AlignedSegment_t
ctypedef utilBMF.HTSUtils.pFastqFile pFastqFile_t
ctypedef utilBMF.HTSUtils.pFastqProxy pFastqProxy_t
ctypedef cython.str cystr

cdef inline cystr GetBS(AlignedSegment_t read):
    return read.opt("BS")

cdef public list int2nuc