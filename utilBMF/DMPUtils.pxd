cimport cython
cimport MawCluster.BCFastq
cimport utilBMF.HTSUtils
from MawCluster.BCFastq cimport cMakeTagComment
ctypedef utilBMF.HTSUtils.pFastqFile pFastqFile_t
ctypedef utilBMF.HTSUtils.pFastqProxy pFastqProxy_t
ctypedef cython.str cystr
ctypedef double double_t
