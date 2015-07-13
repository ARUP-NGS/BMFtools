cimport pysam.calignmentfile
cimport pysam.calignedsegment
cimport cython
from cpython cimport array
ctypedef array.array py_array
ctypedef cython.str cystr

ctypedef pysam.calignedsegment.AlignedSegment AlignedSegment_t
ctypedef pysam.calignmentfile.AlignmentFile AlignmentFile_t

cdef AlignedSegment_t CopyAlignedSegment(AlignedSegment_t template)

cdef class PairwiseAlignmentFile:
    cdef public AlignmentFile_t handle

cdef inline char * PysamToChrInline(char RefID) nogil:
    cdef char * ret
    if RefID == 0:
        ret = '1'
        return ret
    elif RefID == 1:
        ret = '2'
        return ret
    elif RefID == 2:
        ret = '3'
        return ret
    elif RefID == 3:
        ret = '4'
        return ret
    elif RefID == 4:
        ret = '5'
        return ret
    elif RefID == 5:
        ret = '6'
        return ret
    elif RefID == 6:
        ret = '7'
        return ret
    elif RefID == 7:
        ret = '8'
        return ret
    elif RefID == 8:
        ret = '9'
        return ret
    elif RefID == 9:
        ret = '10'
        return ret
    elif RefID == 10:
        ret = '11'
        return ret
    elif RefID == 11:
        ret = '12'
        return ret
    elif RefID == 12:
        ret = '13'
        return ret
    elif RefID == 13:
        ret = '14'
        return ret
    elif RefID == 14:
        ret = '15'
        return ret
    elif RefID == 15:
        ret = '16'
        return ret
    elif RefID == 16:
        ret = '17'
        return ret
    elif RefID == 17:
        ret = '18'
        return ret
    elif RefID == 18:
        ret = '19'
        return ret
    elif RefID == 19:
        ret = '20'
        return ret
    elif RefID == 20:
        ret = '21'
        return ret
    elif RefID == 21:
        ret = '22'
        return ret
    elif RefID == 22:
        ret = 'X'
        return ret
    elif RefID == 23:
        ret = 'Y'
        return ret
    elif RefID == 24:
        ret = 'MT'
        return ret
    elif RefID == 25:
        ret = 'GL000207.1'
        return ret
    elif RefID == 26:
        ret = 'GL000226.1'
        return ret
    elif RefID == 27:
        ret = 'GL000229.1'
        return ret
    elif RefID == 28:
        ret = 'GL000231.1'
        return ret
    elif RefID == 29:
        ret = 'GL000210.1'
        return ret
    elif RefID == 30:
        ret = 'GL000239.1'
        return ret
    elif RefID == 31:
        ret = 'GL000235.1'
        return ret
    elif RefID == 32:
        ret = 'GL000201.1'
        return ret
    elif RefID == 33:
        ret = 'GL000247.1'
        return ret
    elif RefID == 34:
        ret = 'GL000245.1'
        return ret
    elif RefID == 35:
        ret = 'GL000197.1'
        return ret
    elif RefID == 36:
        ret = 'GL000203.1'
        return ret
    elif RefID == 37:
        ret = 'GL000246.1'
        return ret
    elif RefID == 38:
        ret = 'GL000249.1'
        return ret
    elif RefID == 39:
        ret = 'GL000196.1'
        return ret
    elif RefID == 40:
        ret = 'GL000248.1'
        return ret
    elif RefID == 41:
        ret = 'GL000244.1'
        return ret
    elif RefID == 42:
        ret = 'GL000238.1'
        return ret
    elif RefID == 43:
        ret = 'GL000202.1'
        return ret
    elif RefID == 44:
        ret = 'GL000234.1'
        return ret
    elif RefID == 45:
        ret = 'GL000232.1'
        return ret
    elif RefID == 46:
        ret = 'GL000206.1'
        return ret
    elif RefID == 47:
        ret = 'GL000240.1'
        return ret
    elif RefID == 48:
        ret = 'GL000236.1'
        return ret
    elif RefID == 49:
        ret = 'GL000241.1'
        return ret
    elif RefID == 50:
        ret = 'GL000243.1'
        return ret
    elif RefID == 51:
        ret = 'GL000242.1'
        return ret
    elif RefID == 52:
        ret = 'GL000230.1'
        return ret
    elif RefID == 53:
        ret = 'GL000237.1'
        return ret
    elif RefID == 54:
        ret = 'GL000233.1'
        return ret
    elif RefID == 55:
        ret = 'GL000204.1'
        return ret
    elif RefID == 56:
        ret = 'GL000198.1'
        return ret
    elif RefID == 57:
        ret = 'GL000208.1'
        return ret
    elif RefID == 58:
        ret = 'GL000191.1'
        return ret
    elif RefID == 59:
        ret = 'GL000227.1'
        return ret
    elif RefID == 60:
        ret = 'GL000228.1'
        return ret
    elif RefID == 61:
        ret = 'GL000214.1'
        return ret
    elif RefID == 62:
        ret = 'GL000221.1'
        return ret
    elif RefID == 63:
        ret = 'GL000209.1'
        return ret
    elif RefID == 64:
        ret = 'GL000218.1'
        return ret
    elif RefID == 65:
        ret = 'GL000220.1'
        return ret
    elif RefID == 66:
        ret = 'GL000213.1'
        return ret
    elif RefID == 67:
        ret = 'GL000211.1'
        return ret
    elif RefID == 68:
        ret = 'GL000199.1'
        return ret
    elif RefID == 69:
        ret = 'GL000217.1'
        return ret
    elif RefID == 70:
        ret = 'GL000216.1'
        return ret
    elif RefID == 71:
        ret = 'GL000215.1'
        return ret
    elif RefID == 72:
        ret = 'GL000205.1'
        return ret
    elif RefID == 73:
        ret = 'GL000219.1'
        return ret
    elif RefID == 74:
        ret = 'GL000224.1'
        return ret
    elif RefID == 75:
        ret = 'GL000223.1'
        return ret
    elif RefID == 76:
        ret = 'GL000195.1'
        return ret
    elif RefID == 77:
        ret = 'GL000212.1'
        return ret
    elif RefID == 78:
        ret = 'GL000222.1'
        return ret
    elif RefID == 79:
        ret = 'GL000200.1'
        return ret
    elif RefID == 80:
        ret = 'GL000193.1'
        return ret
    elif RefID == 81:
        ret = 'GL000194.1'
        return ret
    elif RefID == 82:
        ret = 'GL000225.1'
        return ret
    else:
        ret = 'GL000192.1'
        return ret
