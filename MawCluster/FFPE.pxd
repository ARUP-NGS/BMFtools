cimport pysam.calignmentfile
cimport pysam.TabProxies
from cpython cimport array
ctypedef array.array py_array
from cython cimport bint
from utilBMF.HTSUtils cimport cystr
ctypedef pysam.calignmentfile.AlignedSegment AlignedSegment_t

cdef AlignedSegment_t AmpliconTrimRead(AlignedSegment_t rec, int primerLen)
cdef dict getFreqDict(pysam.TabProxies.VCFProxy rec)
cdef dict getFreqDictFromInfoDict(dict InfoDict)
cdef dict getCountDictFromInfoDict(dict InfoDict)