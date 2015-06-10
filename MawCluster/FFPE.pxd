cimport pysam.calignmentfile
cimport pysam.TabProxies
from cpython cimport array as c_array
ctypedef c_array.array carray
from cython cimport bint
from utilBMF.HTSUtils cimport cystr
ctypedef pysam.calignmentfile.AlignedSegment cAlignedSegment

cdef cAlignedSegment AmpliconTrimRead(cAlignedSegment rec, int primerLen)
cdef dict getFreqDict(pysam.TabProxies.VCFProxy rec)
cdef dict getFreqDictFromInfoDict(dict InfoDict)
cdef dict getCountDictFromInfoDict(dict InfoDict)