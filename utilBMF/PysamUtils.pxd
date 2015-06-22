cimport pysam.calignmentfile
cimport cython
from cpython cimport array
ctypedef array.array py_array
ctypedef cython.str cystr

ctypedef pysam.calignmentfile.AlignedSegment AlignedSegment_t
ctypedef pysam.calignmentfile.AlignmentFile AlignmentFile_t

cdef AlignedSegment_t CopyAlignedSegment(AlignedSegment_t template)

cdef class PairwiseAlignmentFile:
    cdef public AlignmentFile_t handle

cdef inline cystr PysamToChrInline(char RefID):
    if RefID == 0:
        return '1'
    elif RefID == 1:
        return '2'
    elif RefID == 2:
        return '3'
    elif RefID == 3:
        return '4'
    elif RefID == 4:
        return '5'
    elif RefID == 5:
        return '6'
    elif RefID == 6:
        return '7'
    elif RefID == 7:
        return '8'
    elif RefID == 8:
        return '9'
    elif RefID == 9:
        return '10'
    elif RefID == 10:
        return '11'
    elif RefID == 11:
        return '12'
    elif RefID == 12:
        return '13'
    elif RefID == 13:
        return '14'
    elif RefID == 14:
        return '15'
    elif RefID == 15:
        return '16'
    elif RefID == 16:
        return '17'
    elif RefID == 17:
        return '18'
    elif RefID == 18:
        return '19'
    elif RefID == 19:
        return '20'
    elif RefID == 20:
        return '21'
    elif RefID == 21:
        return '22'
    elif RefID == 22:
        return 'X'
    elif RefID == 23:
        return 'Y'
    elif RefID == 24:
        return 'MT'
    elif RefID == 25:
        return 'GL000207.1'
    elif RefID == 26:
        return 'GL000226.1'
    elif RefID == 27:
        return 'GL000229.1'
    elif RefID == 28:
        return 'GL000231.1'
    elif RefID == 29:
        return 'GL000210.1'
    elif RefID == 30:
        return 'GL000239.1'
    elif RefID == 31:
        return 'GL000235.1'
    elif RefID == 32:
        return 'GL000201.1'
    elif RefID == 33:
        return 'GL000247.1'
    elif RefID == 34:
        return 'GL000245.1'
    elif RefID == 35:
        return 'GL000197.1'
    elif RefID == 36:
        return 'GL000203.1'
    elif RefID == 37:
        return 'GL000246.1'
    elif RefID == 38:
        return 'GL000249.1'
    elif RefID == 39:
        return 'GL000196.1'
    elif RefID == 40:
        return 'GL000248.1'
    elif RefID == 41:
        return 'GL000244.1'
    elif RefID == 42:
        return 'GL000238.1'
    elif RefID == 43:
        return 'GL000202.1'
    elif RefID == 44:
        return 'GL000234.1'
    elif RefID == 45:
        return 'GL000232.1'
    elif RefID == 46:
        return 'GL000206.1'
    elif RefID == 47:
        return 'GL000240.1'
    elif RefID == 48:
        return 'GL000236.1'
    elif RefID == 49:
        return 'GL000241.1'
    elif RefID == 50:
        return 'GL000243.1'
    elif RefID == 51:
        return 'GL000242.1'
    elif RefID == 52:
        return 'GL000230.1'
    elif RefID == 53:
        return 'GL000237.1'
    elif RefID == 54:
        return 'GL000233.1'
    elif RefID == 55:
        return 'GL000204.1'
    elif RefID == 56:
        return 'GL000198.1'
    elif RefID == 57:
        return 'GL000208.1'
    elif RefID == 58:
        return 'GL000191.1'
    elif RefID == 59:
        return 'GL000227.1'
    elif RefID == 60:
        return 'GL000228.1'
    elif RefID == 61:
        return 'GL000214.1'
    elif RefID == 62:
        return 'GL000221.1'
    elif RefID == 63:
        return 'GL000209.1'
    elif RefID == 64:
        return 'GL000218.1'
    elif RefID == 65:
        return 'GL000220.1'
    elif RefID == 66:
        return 'GL000213.1'
    elif RefID == 67:
        return 'GL000211.1'
    elif RefID == 68:
        return 'GL000199.1'
    elif RefID == 69:
        return 'GL000217.1'
    elif RefID == 70:
        return 'GL000216.1'
    elif RefID == 71:
        return 'GL000215.1'
    elif RefID == 72:
        return 'GL000205.1'
    elif RefID == 73:
        return 'GL000219.1'
    elif RefID == 74:
        return 'GL000224.1'
    elif RefID == 75:
        return 'GL000223.1'
    elif RefID == 76:
        return 'GL000195.1'
    elif RefID == 77:
        return 'GL000212.1'
    elif RefID == 78:
        return 'GL000222.1'
    elif RefID == 79:
        return 'GL000200.1'
    elif RefID == 80:
        return 'GL000193.1'
    elif RefID == 81:
        return 'GL000194.1'
    elif RefID == 82:
        return 'GL000225.1'
    else:
        return 'GL000192.1'
