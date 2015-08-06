from pysam.calignmentfile cimport AlignmentFile, IteratorColumnRegion, PileupColumn
from libc.stdint cimport *
from cpython cimport array as c_array
from cython cimport bint
cimport cython

from utilBMF.HTSUtils cimport pPileupRead

# Typedefs
# Cython
ctypedef cython.str cystr
ctypedef c_array.array py_array_t

# Numerical definitions
ctypedef double double_t

# Pysam
ctypedef AlignmentFile AlignmentFile_t
ctypedef PileupColumn PileupColumn_t
ctypedef IteratorColumnRegion IteratorColumnRegion_t

# BMFTools
ctypedef pPileupRead pPileupRead_t


# Class definitions
cdef class SNVAlleleWrangler:
    cdef AlignmentFile_t handle
    cdef public PileupColumn_t column
    cdef public int32_t reference_id, pos, max_insert_size
    cdef public list alleles
    cdef public set insert_sizes
    cdef public dict insert_size_dict
    cdef py_array_t filter_cnf

    cdef inline bint pass_record(self, pPileupRead_t PR, int nuc_index)
    cdef inline void populate(self)
    cdef inline void fast_forward(self)

cdef inline int iabs(int integer) nogil:
    return integer if(integer > 0) else -1 * integer