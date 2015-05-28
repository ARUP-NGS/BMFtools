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

    cdef cython.int getK(self)
    cdef setK(self, cython.int)
    cdef cython.str GetBowtieOutput(self, cython.str)
    cdef FillMap(self, list)
    cdef list GetUniqueKmers(self, list)
    cdef list GetKmerList(self, list)