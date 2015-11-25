# cython: boundscheck=False, wraparound=False
# cython: cdivision=True, initializedcheck=False
"""
A module composed solely of precomputed inline functions. This way, when the
modules that these originally came from are modified and rebuilt, at least
the time to compile these portions can be saved.
"""
cimport cython
ctypedef cython.str cystr
cimport numpy as np
ctypedef np.int32_t np_int32_t

from libc.string cimport strcmp
from libc.stdint cimport int8_t, int16_t
from libc.stdio cimport sprintf


cdef inline char * opLenToStr(char op, int opLen) nogil:
    cdef char[5] ret
    if(op == 68):
        sprintf("%iD", ret, opLen)
    elif(op == 77):
        sprintf("%iM", ret, opLen)
    elif(op == 73):
        sprintf("%iI", ret, opLen)
    elif(op == 83):
        sprintf("%iS", ret, opLen)
    else:
        sprintf("%iN", ret, opLen)
    return ret


cdef inline char * CigarStrInline(char character) nogil:
    cdef char * ret
    if character == 77:
        ret = "M"
    elif character == 73:
        ret = "I"
    elif character == 68:
        ret = "D"
    elif character == 78:
        ret = "N"
    elif character == 83:
        ret = "S"
    elif character == 72:
        ret = "H"
    elif character == 80:
        ret = "P"
    elif character == 61:
        ret = "="
    else:
        ret = "X"
    return ret


cdef inline char * Num2Nuc(int number) nogil:
    cdef char * ret
    if(number == 0):
        ret = "A"
    elif(number == 1):
        ret = "C"
    elif(number == 2):
        ret = "G"
    elif(number == 3):
        ret = "T"
    else:
        ret = "N"
    return ret


cdef inline char Num2NucN(int number) nogil:
    if number == 0:
        return 65
    elif number == 1:
        return 67
    elif number == 2:
        return 71
    elif number == 3:
        return 84
    else:
        return 78


cdef inline char Nuc2NumN(char character) nogil:
    if(character == 84):
        return 3
    elif(character == 71):
        return 2
    elif(character == 67):
        return 1
    elif character == 65:
        return 0
    else:
        return 4


cdef inline char Nuc2Num(char character) nogil:
    if(character == 84):
        return 3
    elif(character == 71):
        return 2
    elif(character == 67):
        return 1
    else:
        return 0


cdef inline char * chrACGNTInline(char character) nogil:
    cdef char * ret
    if(character == 65):
        ret = "A"
    elif(character == 67):
        ret = "C"
    elif(character == 71):
        ret = "G"
    elif(character == 84):
        ret = "T"
    elif(character == 78):
        ret = "N"
    else:
        ret = "X"
    return ret


cdef inline char * chrInline(char character) nogil:
    cdef char * ret
    if(character == 33):
            ret = '!'
    elif(character == 34):
            ret = '"'
    elif(character == 35):
            ret = '#'
    elif(character == 36):
            ret = '$'
    elif(character == 37):
            ret = '%'
    elif(character == 38):
            ret = '&'
    elif(character == 39):
            ret = '\''
    elif(character == 40):
            ret = '('
    elif(character == 41):
            ret = ')'
    elif(character == 42):
            ret = '*'
    elif(character == 43):
            ret = '+'
    elif(character == 44):
            ret = ','
    elif(character == 45):
            ret = '-'
    elif(character == 46):
            ret = '.'
    elif(character == 47):
            ret = '/'
    elif(character == 48):
            ret = '0'
    elif(character == 49):
            ret = '1'
    elif(character == 50):
            ret = '2'
    elif(character == 51):
            ret = '3'
    elif(character == 52):
            ret = '4'
    elif(character == 53):
            ret = '5'
    elif(character == 54):
            ret = '6'
    elif(character == 55):
            ret = '7'
    elif(character == 56):
            ret = '8'
    elif(character == 57):
            ret = '9'
    elif(character == 58):
            ret = ':'
    elif(character == 59):
            ret = ';'
    elif(character == 60):
            ret = '<'
    elif(character == 61):
            ret = '='
    elif(character == 62):
            ret = '>'
    elif(character == 63):
            ret = '?'
    elif(character == 64):
            ret = '@'
    elif(character == 65):
            ret = 'A'
    elif(character == 66):
            ret = 'B'
    elif(character == 67):
            ret = 'C'
    elif(character == 68):
            ret = 'D'
    elif(character == 69):
            ret = 'E'
    elif(character == 70):
            ret = 'F'
    elif(character == 71):
            ret = 'G'
    elif(character == 72):
            ret = 'H'
    elif(character == 73):
            ret = 'I'
    elif(character == 74):
            ret = 'J'
    elif(character == 75):
            ret = 'K'
    elif(character == 76):
            ret = 'L'
    elif(character == 77):
            ret = 'M'
    elif(character == 78):
            ret = 'N'
    elif(character == 79):
            ret = 'O'
    elif(character == 80):
            ret = 'P'
    elif(character == 81):
            ret = 'Q'
    elif(character == 82):
            ret = 'R'
    elif(character == 83):
            ret = 'S'
    elif(character == 84):
            ret = 'T'
    elif(character == 85):
            ret = 'U'
    elif(character == 86):
            ret = 'V'
    elif(character == 87):
            ret = 'W'
    elif(character == 88):
            ret = 'X'
    elif(character == 89):
            ret = 'Y'
    elif(character == 90):
            ret = 'Z'
    elif(character == 91):
            ret = '['
    elif(character == 92):
            ret = '\\'
    elif(character == 93):
            ret = ']'
    elif(character == 94):
            ret = '^'
    elif(character == 95):
            ret = '_'
    elif(character == 96):
            ret = '`'
    elif(character == 97):
            ret = 'a'
    elif(character == 98):
            ret = 'b'
    elif(character == 99):
            ret = 'c'
    elif(character == 100):
            ret = 'd'
    elif(character == 101):
            ret = 'e'
    elif(character == 102):
            ret = 'f'
    elif(character == 103):
            ret = 'g'
    elif(character == 104):
            ret = 'h'
    elif(character == 105):
            ret = 'i'
    elif(character == 106):
            ret = 'j'
    elif(character == 107):
            ret = 'k'
    elif(character == 108):
            ret = 'l'
    elif(character == 109):
            ret = 'm'
    elif(character == 110):
            ret = 'n'
    elif(character == 111):
            ret = 'o'
    elif(character == 112):
            ret = 'p'
    elif(character == 113):
            ret = 'q'
    elif(character == 114):
            ret = 'r'
    elif(character == 115):
            ret = 's'
    elif(character == 116):
            ret = 't'
    elif(character == 117):
            ret = 'u'
    elif(character == 118):
            ret = 'v'
    elif(character == 119):
            ret = 'w'
    elif(character == 120):
            ret = 'x'
    elif(character == 121):
            ret = 'y'
    elif(character == 122):
            ret = 'z'
    elif(character == 123):
            ret = '{'
    elif(character == 124):
            ret = '|'
    elif(character == 125):
            ret = '}'
    else:
            ret = '~'
    return ret


cdef inline char CigarOpToCigarChar(char character) nogil:
    if(character == 0):
        return 77
    elif(character == 1):
        return 73
    elif(character == 2):
        return 68
    elif(character == 3):
        return 78
    elif(character == 4):
        return 83
    elif(character == 5):
        return 72
    elif(character == 6):
        return 80
    elif(character == 7):
        return 61
    elif(character == 8):
        return 88
    elif(character == 73):
        return 1
    elif(character == 77):
        return 0
    elif(character == 78):
        return 3
    elif(character == 80):
        return 6
    elif(character == 72):
        return 5
    elif(character == 83):
        return 4
    elif(character == 88):
        return 8
    elif(character == 68):
        return 2
    else:
        return 7


cdef inline char ChrToRefIDInline(char * contig) nogil:
    if(strcmp(contig, '1') == 0):
        return 0
    elif(strcmp(contig, '2') == 0):
        return 1
    elif(strcmp(contig, '3') == 0):
        return 2
    elif(strcmp(contig, '4') == 0):
        return 3
    elif(strcmp(contig, '5') == 0):
        return 4
    elif(strcmp(contig, '6') == 0):
        return 5
    elif(strcmp(contig, '7') == 0):
        return 6
    elif(strcmp(contig, '8') == 0):
        return 7
    elif(strcmp(contig, '9') == 0):
        return 8
    elif(strcmp(contig, '10') == 0):
        return 9
    elif(strcmp(contig, '11') == 0):
        return 10
    elif(strcmp(contig, '12') == 0):
        return 11
    elif(strcmp(contig, '13') == 0):
        return 12
    elif(strcmp(contig, '14') == 0):
        return 13
    elif(strcmp(contig, '15') == 0):
        return 14
    elif(strcmp(contig, '16') == 0):
        return 15
    elif(strcmp(contig, '17') == 0):
        return 16
    elif(strcmp(contig, '18') == 0):
        return 17
    elif(strcmp(contig, '19') == 0):
        return 18
    elif(strcmp(contig, '20') == 0):
        return 19
    elif(strcmp(contig, '21') == 0):
        return 20
    elif(strcmp(contig, '22') == 0):
        return 21
    elif(strcmp(contig, 'X') == 0):
        return 22
    elif(strcmp(contig, 'Y') == 0):
        return 23
    elif(strcmp(contig, 'MT') == 0):
        return 24
    elif(strcmp(contig, 'GL000207.1') == 0):
        return 25
    elif(strcmp(contig, 'GL000226.1') == 0):
        return 26
    elif(strcmp(contig, 'GL000229.1') == 0):
        return 27
    elif(strcmp(contig, 'GL000231.1') == 0):
        return 28
    elif(strcmp(contig, 'GL000210.1') == 0):
        return 29
    elif(strcmp(contig, 'GL000239.1') == 0):
        return 30
    elif(strcmp(contig, 'GL000235.1') == 0):
        return 31
    elif(strcmp(contig, 'GL000201.1') == 0):
        return 32
    elif(strcmp(contig, 'GL000247.1') == 0):
        return 33
    elif(strcmp(contig, 'GL000245.1') == 0):
        return 34
    elif(strcmp(contig, 'GL000197.1') == 0):
        return 35
    elif(strcmp(contig, 'GL000203.1') == 0):
        return 36
    elif(strcmp(contig, 'GL000246.1') == 0):
        return 37
    elif(strcmp(contig, 'GL000249.1') == 0):
        return 38
    elif(strcmp(contig, 'GL000196.1') == 0):
        return 39
    elif(strcmp(contig, 'GL000248.1') == 0):
        return 40
    elif(strcmp(contig, 'GL000244.1') == 0):
        return 41
    elif(strcmp(contig, 'GL000238.1') == 0):
        return 42
    elif(strcmp(contig, 'GL000202.1') == 0):
        return 43
    elif(strcmp(contig, 'GL000234.1') == 0):
        return 44
    elif(strcmp(contig, 'GL000232.1') == 0):
        return 45
    elif(strcmp(contig, 'GL000206.1') == 0):
        return 46
    elif(strcmp(contig, 'GL000240.1') == 0):
        return 47
    elif(strcmp(contig, 'GL000236.1') == 0):
        return 48
    elif(strcmp(contig, 'GL000241.1') == 0):
        return 49
    elif(strcmp(contig, 'GL000243.1') == 0):
        return 50
    elif(strcmp(contig, 'GL000242.1') == 0):
        return 51
    elif(strcmp(contig, 'GL000230.1') == 0):
        return 52
    elif(strcmp(contig, 'GL000237.1') == 0):
        return 53
    elif(strcmp(contig, 'GL000233.1') == 0):
        return 54
    elif(strcmp(contig, 'GL000204.1') == 0):
        return 55
    elif(strcmp(contig, 'GL000198.1') == 0):
        return 56
    elif(strcmp(contig, 'GL000208.1') == 0):
        return 57
    elif(strcmp(contig, 'GL000191.1') == 0):
        return 58
    elif(strcmp(contig, 'GL000227.1') == 0):
        return 59
    elif(strcmp(contig, 'GL000228.1') == 0):
        return 60
    elif(strcmp(contig, 'GL000214.1') == 0):
        return 61
    elif(strcmp(contig, 'GL000221.1') == 0):
        return 62
    elif(strcmp(contig, 'GL000209.1') == 0):
        return 63
    elif(strcmp(contig, 'GL000218.1') == 0):
        return 64
    elif(strcmp(contig, 'GL000220.1') == 0):
        return 65
    elif(strcmp(contig, 'GL000213.1') == 0):
        return 66
    elif(strcmp(contig, 'GL000211.1') == 0):
        return 67
    elif(strcmp(contig, 'GL000199.1') == 0):
        return 68
    elif(strcmp(contig, 'GL000217.1') == 0):
        return 69
    elif(strcmp(contig, 'GL000216.1') == 0):
        return 70
    elif(strcmp(contig, 'GL000215.1') == 0):
        return 71
    elif(strcmp(contig, 'GL000205.1') == 0):
        return 72
    elif(strcmp(contig, 'GL000219.1') == 0):
        return 73
    elif(strcmp(contig, 'GL000224.1') == 0):
        return 74
    elif(strcmp(contig, 'GL000223.1') == 0):
        return 75
    elif(strcmp(contig, 'GL000195.1') == 0):
        return 76
    elif(strcmp(contig, 'GL000212.1') == 0):
        return 77
    elif(strcmp(contig, 'GL000222.1') == 0):
        return 78
    elif(strcmp(contig, 'GL000200.1') == 0):
        return 79
    elif(strcmp(contig, 'GL000193.1') == 0):
        return 80
    elif(strcmp(contig, 'GL000194.1') == 0):
        return 81
    elif(strcmp(contig, 'GL000225.1') == 0):
        return 82
    else:
        return 83


cdef inline int16_t CONTEXT_TO_ARRAY_POS(char * context) nogil:
    if(strcmp(context, 'AA') == 0):
        return 0
    elif(strcmp(context, 'AC') == 0):
        return 1
    elif(strcmp(context, 'AG') == 0):
        return 2
    elif(strcmp(context, 'AT') == 0):
        return 3
    elif(strcmp(context, 'CA') == 0):
        return 4
    elif(strcmp(context, 'CC') == 0):
        return 5
    elif(strcmp(context, 'CG') == 0):
        return 6
    elif(strcmp(context, 'CT') == 0):
        return 7
    elif(strcmp(context, 'GA') == 0):
        return 8
    elif(strcmp(context, 'GC') == 0):
        return 9
    elif(strcmp(context, 'GG') == 0):
        return 10
    elif(strcmp(context, 'GT') == 0):
        return 11
    elif(strcmp(context, 'TA') == 0):
        return 12
    elif(strcmp(context, 'TC') == 0):
        return 13
    elif(strcmp(context, 'TG') == 0):
        return 14
    elif(strcmp(context, 'TT') == 0):
        return 15
    else:
        return -1
