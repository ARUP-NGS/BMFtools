cimport cython
cimport numpy as np
cimport pysam.calignmentfile
cimport pysam.cfaidx
cimport utilBMF.HTSUtils
ctypedef utilBMF.HTSUtils.pFastqFile pFastqFile_t
ctypedef utilBMF.HTSUtils.pFastqProxy pFastqProxy_t
from utilBMF.HTSUtils cimport cystr, chr2ph, chr2phStr, int2Str, ph2chrDict
from numpy cimport ndarray

cdef cystr compareFqRecsFast(list R,
                                  cystr name=?,
                                  int famLimit=?,
                                  dict chr2ph=?,
                                  dict letterNumDict=?,
                                  dict ph2chrDict=?,
                                  object ccopy=?,
                                  object npchararray=?,
                                  object oagqual=?,
                                  object oagseq=?,
                                  object partialnpchar=?,
                                  object ph2chr=?,
                                  object chr2phStr=?,
                                  object int2Str=?)
cdef cystr compareFqRecsFqPrx(list R,
                                   cystr name=?,
                                   float stringency=?,
                                   int famLimit=?,
                                   cython.bint keepFails=?,
                                   object oagseq=?,
                                   dict chr2ph=?,
                                   dict ph2chrDict=?,
                                   object ph2chr=?,
                                   object int2Str=?)
cpdef cystr cFRF_helper(list R, cystr name=?)
cpdef cystr cFRP_helper(list R, cystr name=?)
cpdef cystr MakeTagComment(cystr saltedBS, pFastqProxy_t rec, int)
cdef cystr cMakeTagComment(cystr saltedBS, pFastqProxy_t rec, int)

cdef public dict letterNumDict