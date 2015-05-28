cimport cython

cdef class KmerFetcher(object):
    """
    Contains parameters for finding representative kmers.
    I want to permit the mixing of kmers of different sizes, so I have
    set a baseK, which is the K that we start using.
    """
    cdef cython.int k

    cdef public cython.str ref
    cdef public cython.int mismatches, minMQ, padding
    cdef public dict HashMap

    cdef public cython.int getK(self)
    cdef public setK(self, cython.int)
    cpdef cython.str GetBowtieOutput(self, cython.str)
    cpdef public FillMap(self, list)
    cpdef public list GetUniqueKmers(self, list)
    cdef public list GetKmerList(self, list)