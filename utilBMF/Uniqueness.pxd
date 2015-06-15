cimport cython
from cython cimport bint
from utilBMF.HTSUtils cimport cystr

cdef class KmerFetcher(object):
    """
    Contains parameters for finding representative kmers.
    I want to permit the mixing of kmers of different sizes, so I have
    set a baseK, which is the K that we start using.
    """
    cdef readonly cython.int k

    cdef readonly cystr ref
    cdef public cython.int mismatches, minMQ, padding
    cdef public dict HashMap, FullMap

    cdef public cython.int getK(self)
    cdef public setK(self, cython.int)
    cpdef cystr getFastqString(self, list)
    cpdef cystr getOutputString(self, list, str aligner=?)
    cpdef public FillMap(self, list)
    cpdef public list GetUniqueKmers(self, list)
    cpdef FMfrombed(self, cystr)


cdef class RefKmer(object):
    cdef readonly cystr contig, seq
    cdef readonly cython.int len, pos

cpdef cystr SequenceToFakeFq(cystr seq)
cpdef list GetKmersToCheck(cystr ref, int k=?, list bedline=?,
                           int padding=?)
cpdef cystr FastqStrFromKmerList(list kmerList)
cpdef list GetRepKmersBwt(cystr ref,
                          int k=?,
                          list bedline=?,
                          int padding=?,
                          int seedlen=?,
                          int mismatches=?,
                          int minMQ=?,
                          bint useBowtie=?)
cpdef cystr BowtieFqToStr(cystr fqStr, cystr ref=?,
                          int seed=?, int mismatches=?)

cdef dict mmDict
