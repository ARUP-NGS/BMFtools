cimport cython
cimport numpy as np
cimport pysam.calignmentfile
cimport pysam.cfaidx
cimport utilBMF.HTSUtils
ctypedef utilBMF.HTSUtils.pFastqFile pFastqFile_t
ctypedef utilBMF.HTSUtils.pFastqProxy pFastqProxy_t
from utilBMF.HTSUtils cimport chr2ph, chr2phStr, int2Str, ph2chrDict
from utilBMF.cstring cimport cs_to_ph
from utilBMF.Inliners cimport Num2Nuc, Nuc2Num, ph2chrInline, int2strInline, chr2phInline, chr2phImplicit
from numpy cimport ndarray
from cpython cimport array
import cython.str
ctypedef cython.str cystr
ctypedef array.array py_array
ctypedef np.int32_t np_int32_t

cdef cystr cCompareFqRecsFast(list R,
                              cystr name=?,
                              int famLimit=?)
cdef cystr compareFqRecsFqPrx(list R,
                              cystr name=?,
                              float stringency=?,
                              int famLimit=?,
                              cython.bint keepFails=?,
                              object oagseq=?)
cpdef cystr pCompareFqRecsFast(list R, cystr name=?)
cpdef cystr cFRP_helper(list R, cystr name=?)
cpdef cystr MakeTagComment(cystr saltedBS, pFastqProxy_t rec, int)
cdef cystr cMakeTagComment(cystr saltedBS, pFastqProxy_t rec, int)

cdef cystr cQualArr2QualStr(ndarray[np_int32_t, ndim=1] qualArr)
cpdef cystr QualArr2QualStr(ndarray[np_int32_t, ndim=1] qualArr)

cpdef cystr QualArr2PVString(ndarray[np_int32_t, ndim=1] qualArr)
cdef cystr cQualArr2PVString(ndarray[np_int32_t, ndim=1] qualArr)

