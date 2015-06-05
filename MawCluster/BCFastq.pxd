cimport cython
cimport numpy as np
cimport pysam.calignmentfile
cimport pysam.cfaidx
cimport utilBMF.HTSUtils
ctypedef utilBMF.HTSUtils.pFastqFile pFastqFile_t
ctypedef utilBMF.HTSUtils.pFastqProxy pFastqProxy_t
from numpy cimport ndarray

cdef cython.str compareFqRecsFast(list R,
                                  cython.str name=?,
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
cdef cython.str compareFqRecsFqPrx(list R,
                                   cython.str name=?,
                                   float stringency=?,
                                   int famLimit=?,
                                   cython.bint keepFails=?,
                                   object oagseq=?,
                                   dict chr2ph=?,
                                   dict ph2chrDict=?,
                                   object ph2chr=?,
                                   object int2Str=?)
cpdef cython.str cFRF_helper(list R, cython.str name=?)
cpdef cython.str cFRP_helper(list R, cython.str name=?)
