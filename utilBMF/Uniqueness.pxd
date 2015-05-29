cimport cython

cdef class KmerFetcher(object):
    """
    Contains parameters for finding representative kmers.
    I want to permit the mixing of kmers of different sizes, so I have
    set a baseK, which is the K that we start using.
    """
    cdef readonly cython.int k

    cdef readonly cython.str ref
    cdef public cython.int mismatches, minMQ, padding
    cdef public dict HashMap

    cdef public cython.int getK(self)
    cdef public setK(self, cython.int)
    cpdef cython.str getFastqString(self, list)
    cpdef cython.str getOutputString(self, list)
    cpdef public FillMap(self, list)
    cpdef public list GetUniqueKmers(self, list)
    cpdef FMfrombed(self, cython.str)


cdef class RefKmer(object):
    cdef readonly cython.str contig, seq
    cdef readonly cython.int len, pos