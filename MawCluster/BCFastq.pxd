cimport cython
cimport utilBMF.HTSUtils
# from utilBMF._re2 cimport cpp_string
ctypedef utilBMF.HTSUtils.pFastqProxy pFq
ctypedef utilBMF.HTSUtils.pFastqFile pFastqFile_t

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
                                  object fromQS=?)
cdef cython.str compareFqRecsFqPrx(list R,
                                   cython.str name=?,
                                   float stringency=?,
                                   int famLimit=?, cython.bint keepFails=?,
                                   object oagseq=?,
                                   dict chr2ph=?,
                                   dict ph2chrDict=?,
                                   object ph2chr=?,
                                   object fromQS=?)
cpdef cython.str cFRF_helper(list R, cython.str name=?)
cpdef cython.str cFRP_helper(list R, cython.str name=?)
