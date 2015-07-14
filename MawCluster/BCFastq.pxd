cimport cython
cimport numpy as np
cimport pysam.calignmentfile
cimport pysam.cfaidx
cimport utilBMF.HTSUtils
from cpython cimport array as c_array
from libc.stdlib cimport malloc, free, realloc
from libc.stdint cimport int8_t
from libc.string cimport memcpy
from numpy cimport ndarray, uint8_t
from utilBMF.cstring cimport cs_to_ph, cs_to_ia, PH2CHR_TRANS
from utilBMF.HTSUtils cimport chr2ph, chr2phStr, int2Str, ph2chrDict
from utilBMF.Inliners cimport Num2Nuc, Nuc2Num
ctypedef c_array.array py_array
ctypedef cython.str cystr
ctypedef SeqQual SeqQual_t
ctypedef utilBMF.HTSUtils.pFastqFile pFastqFile_t
ctypedef utilBMF.HTSUtils.pFastqProxy pFastqProxy_t
ctypedef int int32_t


# METHODS
cdef cystr cCompareFqRecsFast(list R, cystr name=?, double minPVFrac=?,
                              double minFAFrac=?, double minMaxFA=?)
cpdef cystr pCompareFqRecsFast(list R, cystr name=?)
cpdef cystr MakeTagComment(cystr saltedBS, pFastqProxy_t rec, int)

cdef cystr cQualArr2QualStr(ndarray[int32_t, ndim=1] qualArr)
cpdef cystr QualArr2QualStr(ndarray[int32_t, ndim=1] qualArr)

cpdef cystr QualArr2PVString(ndarray[int32_t, ndim=1] qualArr)
cdef cystr cQualArr2PVString(ndarray[int32_t, ndim=1] qualArr)
cdef cystr cQualArr2FAString(ndarray[int32_t, ndim=1] qualArr)


cdef inline bint BarcodePasses(cystr barcode, int hpLimit)


cdef inline cystr cMakeTagComment(cystr saltedBS,
                                  pFastqProxy_t rec, int hpLimit):
    cdef bint PASS = BarcodePasses(saltedBS, hpLimit)
    if(PASS):
        return "~#!#~" + rec.comment + "|FP=1|BS=" + saltedBS
    else:
        return "~#!#~" + rec.comment + "|FP=0|BS=" + saltedBS


# CONSTANTS
cdef public dict Num2NucDict
cdef public cystr ARGMAX_TRANSLATE_STRING
cdef public py_array nucs
cdef char size32 = 4
cdef char sizedouble = 8

# STRUCTS
cdef struct SeqQual:
    np.int8_t * Seq
    int32_t * Agree
    int32_t * Qual
    size_t length
