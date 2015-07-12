cimport pysam.calignmentfile
cimport pysam.calignedsegment
cimport pysam.ctabixproxies
from cpython cimport array
ctypedef array.array py_array
from cython cimport bint
from utilBMF.HTSUtils cimport cystr
ctypedef pysam.calignedsegment.AlignedSegment AlignedSegment_t

cdef AlignedSegment_t AmpliconTrimRead(AlignedSegment_t rec, int primerLen)
cdef dict getFreqDict(pysam.ctabixproxies.VCFProxy rec)
cdef dict getFreqDictFromInfoDict(dict InfoDict)
cdef dict getCountDictFromInfoDict(dict InfoDict)
