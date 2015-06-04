cimport pysam.cfaidx
cimport numpy as np
cimport utilBMF.HTSUtils
ctypedef utilBMF.HTSUtils.pFastqFile pFastqFile_t
ctypedef utilBMF.HTSUtils.pFastqProxy pFastqProxy_t
cdef np.ndarray[double] GetFamSizeStats_(pFastqFile_t FqHandle)
