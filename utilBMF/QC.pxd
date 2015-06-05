cimport pysam.cfaidx
cimport numpy as np
from numpy cimport ndarray
cimport utilBMF.HTSUtils
ctypedef utilBMF.HTSUtils.pFastqFile pFastqFile_t
ctypedef utilBMF.HTSUtils.pFastqProxy pFastqProxy_t
cdef ndarray[double] GetFamSizeStats_(pFastqFile_t FqHandle)
cdef ndarray[np.float64_t, ndim = 1] InsertSizeArray_(pysam.calignmentfile.AlignmentFile)