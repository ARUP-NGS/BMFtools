cimport pysam.calignmentfile
cimport pysam.TabProxies
from cpython cimport array as c_array
ctypedef c_array.array carray

cdef dict getFreqDict(pysam.TabProxies.VCFProxy rec)
cdef dict getFreqDictFromInfoDict(dict InfoDict)