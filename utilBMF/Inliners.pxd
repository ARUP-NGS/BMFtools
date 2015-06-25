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


cdef inline cystr opLenToStr(char op, int opLen):
    if(op == 68):
        return "%sD" % opLen
    elif(op == 77):
        return "%sM" % opLen
    elif(op == 73):
        return "%sI" % opLen
    elif(op == 83):
        return "%sS" % opLen
    else:
        return ""


cdef inline cystr CigarStrInline(char character):
    if character == 77:
        return "M"
    elif character == 73:
        return "I"
    elif character == 68:
        return "D"
    elif character == 78:
        return "N"
    elif character == 83:
        return "S"
    elif character == 72:
        return "H"
    elif character == 80:
        return "P"
    elif character == 61:
        return "="
    else:
        return "X"


cdef inline cystr RevCmpChar(cystr character):
    if(character == "A"):
        return "T"
    elif(character == "C"):
        return "G"
    elif(character == "G"):
        return "C"
    elif(character == "T"):
        return "A"
    else:
        return "N"


cdef inline cystr RevCmpInt(char character):
    if(character == 65):
        return "T"
    elif(character == 67):
        return "G"
    elif(character == 84):
        return "A"
    elif(character == 71):
        return "C"
    else:
        return "N"


cdef inline char RevCmpToChar(char character) nogil:
    if(character == 65):
        return 84
    elif(character == 67):
        return 71
    elif(character == 84):
        return 65
    elif(character == 71):
        return 57
    else:
        return 78


cdef inline cystr Num2Nuc(int number):
    if(number == 0):
        return "A"
    elif(number == 1):
        return "C"
    elif(number == 2):
        return "G"
    elif(number == 3):
        return "T"
    else:
        return "N"


cdef inline char Nuc2Num(char character) nogil:
    if(character == 84):
        return 3
    elif(character == 71):
        return 2
    elif(character == 67):
        return 1
    else:
        return 0


cdef inline cystr ph2chrInline(np_int32_t phInput):
    if(phInput == 0):
        return '!'
    elif(phInput == 1):
        return '"'
    elif(phInput == 2):
        return '#'
    elif(phInput == 3):
        return '$'
    elif(phInput == 4):
        return '%'
    elif(phInput == 5):
        return '&'
    elif(phInput == 6):
        return '\''
    elif(phInput == 7):
        return '('
    elif(phInput == 8):
        return ')'
    elif(phInput == 9):
        return '*'
    elif(phInput == 10):
        return '+'
    elif(phInput == 11):
        return ','
    elif(phInput == 12):
        return '-'
    elif(phInput == 13):
        return '.'
    elif(phInput == 14):
        return '/'
    elif(phInput == 15):
        return '0'
    elif(phInput == 16):
        return '1'
    elif(phInput == 17):
        return '2'
    elif(phInput == 18):
        return '3'
    elif(phInput == 19):
        return '4'
    elif(phInput == 20):
        return '5'
    elif(phInput == 21):
        return '6'
    elif(phInput == 22):
        return '7'
    elif(phInput == 23):
        return '8'
    elif(phInput == 24):
        return '9'
    elif(phInput == 25):
        return ':'
    elif(phInput == 26):
        return ';'
    elif(phInput == 27):
        return '<'
    elif(phInput == 28):
        return '='
    elif(phInput == 29):
        return '>'
    elif(phInput == 30):
        return '?'
    elif(phInput == 31):
        return '@'
    elif(phInput == 32):
        return 'A'
    elif(phInput == 33):
        return 'B'
    elif(phInput == 34):
        return 'C'
    elif(phInput == 35):
        return 'D'
    elif(phInput == 36):
        return 'E'
    elif(phInput == 37):
        return 'F'
    elif(phInput == 38):
        return 'G'
    elif(phInput == 39):
        return 'H'
    elif(phInput == 40):
        return 'I'
    elif(phInput == 41):
        return 'J'
    elif(phInput == 42):
        return 'K'
    elif(phInput == 43):
        return 'L'
    elif(phInput == 44):
        return 'M'
    elif(phInput == 45):
        return 'N'
    elif(phInput == 46):
        return 'O'
    elif(phInput == 47):
        return 'P'
    elif(phInput == 48):
        return 'Q'
    elif(phInput == 49):
        return 'R'
    elif(phInput == 50):
        return 'S'
    elif(phInput == 51):
        return 'T'
    elif(phInput == 52):
        return 'U'
    elif(phInput == 53):
        return 'V'
    elif(phInput == 54):
        return 'W'
    elif(phInput == 55):
        return 'X'
    elif(phInput == 56):
        return 'Y'
    elif(phInput == 57):
        return 'Z'
    elif(phInput == 58):
        return '['
    elif(phInput == 59):
        return '\\'
    elif(phInput == 60):
        return ']'
    elif(phInput == 61):
        return '^'
    elif(phInput == 62):
        return '_'
    elif(phInput == 63):
        return '`'
    elif(phInput == 64):
        return 'a'
    elif(phInput == 65):
        return 'b'
    elif(phInput == 66):
        return 'c'
    elif(phInput == 67):
        return 'd'
    elif(phInput == 68):
        return 'e'
    elif(phInput == 69):
        return 'f'
    elif(phInput == 70):
        return 'g'
    elif(phInput == 71):
        return 'h'
    elif(phInput == 72):
        return 'i'
    elif(phInput == 73):
        return 'j'
    elif(phInput == 74):
        return 'k'
    elif(phInput == 75):
        return 'l'
    elif(phInput == 76):
        return 'm'
    elif(phInput == 77):
        return 'n'
    elif(phInput == 78):
        return 'o'
    elif(phInput == 79):
        return 'p'
    elif(phInput == 80):
        return 'q'
    elif(phInput == 81):
        return 'r'
    elif(phInput == 82):
        return 's'
    elif(phInput == 83):
        return 't'
    elif(phInput == 84):
        return 'u'
    elif(phInput == 85):
        return 'v'
    elif(phInput == 86):
        return 'w'
    elif(phInput == 87):
        return 'x'
    elif(phInput == 88):
        return 'y'
    elif(phInput == 89):
        return 'z'
    elif(phInput == 90):
        return '{'
    elif(phInput == 91):
        return '|'
    elif(phInput == 92):
        return '}'
    else:
        return '~'


cdef inline char chr2phInline(cystr character):
    if(character == '$'):
        return 3
    elif(character == '('):
        return 7
    elif(character == ','):
        return 11
    elif(character == '0'):
        return 15
    elif(character == '4'):
        return 19
    elif(character == '8'):
        return 23
    elif(character == '<'):
        return 27
    elif(character == '@'):
        return 31
    elif(character == 'D'):
        return 35
    elif(character == 'H'):
        return 39
    elif(character == 'L'):
        return 43
    elif(character == 'P'):
        return 47
    elif(character == 'T'):
        return 51
    elif(character == 'X'):
        return 55
    elif(character == '\\'):
        return 59
    elif(character == '`'):
        return 63
    elif(character == 'd'):
        return 67
    elif(character == 'h'):
        return 71
    elif(character == 'l'):
        return 75
    elif(character == 'p'):
        return 79
    elif(character == 't'):
        return 83
    elif(character == 'x'):
        return 87
    elif(character == '|'):
        return 91
    elif(character == '#'):
        return 2
    elif(character == '\''):
        return 6
    elif(character == '+'):
        return 10
    elif(character == '/'):
        return 14
    elif(character == '3'):
        return 18
    elif(character == '7'):
        return 22
    elif(character == ';'):
        return 26
    elif(character == '?'):
        return 30
    elif(character == 'C'):
        return 34
    elif(character == 'G'):
        return 38
    elif(character == 'K'):
        return 42
    elif(character == 'O'):
        return 46
    elif(character == 'S'):
        return 50
    elif(character == 'W'):
        return 54
    elif(character == '['):
        return 58
    elif(character == '_'):
        return 62
    elif(character == 'c'):
        return 66
    elif(character == 'g'):
        return 70
    elif(character == 'k'):
        return 74
    elif(character == 'o'):
        return 78
    elif(character == 's'):
        return 82
    elif(character == 'w'):
        return 86
    elif(character == '{'):
        return 90
    elif(character == '"'):
        return 1
    elif(character == '&'):
        return 5
    elif(character == '*'):
        return 9
    elif(character == '.'):
        return 13
    elif(character == '2'):
        return 17
    elif(character == '6'):
        return 21
    elif(character == ':'):
        return 25
    elif(character == '>'):
        return 29
    elif(character == 'B'):
        return 33
    elif(character == 'F'):
        return 37
    elif(character == 'J'):
        return 41
    elif(character == 'N'):
        return 45
    elif(character == 'R'):
        return 49
    elif(character == 'V'):
        return 53
    elif(character == 'Z'):
        return 57
    elif(character == '^'):
        return 61
    elif(character == 'b'):
        return 65
    elif(character == 'f'):
        return 69
    elif(character == 'j'):
        return 73
    elif(character == 'n'):
        return 77
    elif(character == 'r'):
        return 81
    elif(character == 'v'):
        return 85
    elif(character == 'z'):
        return 89
    elif(character == '!'):
        return 0
    elif(character == '%'):
        return 4
    elif(character == ')'):
        return 8
    elif(character == '-'):
        return 12
    elif(character == '1'):
        return 16
    elif(character == '5'):
        return 20
    elif(character == '9'):
        return 24
    elif(character == '='):
        return 28
    elif(character == 'A'):
        return 32
    elif(character == 'E'):
        return 36
    elif(character == 'I'):
        return 40
    elif(character == 'M'):
        return 44
    elif(character == 'Q'):
        return 48
    elif(character == 'U'):
        return 52
    elif(character == 'Y'):
        return 56
    elif(character == ']'):
        return 60
    elif(character == 'a'):
        return 64
    elif(character == 'e'):
        return 68
    elif(character == 'i'):
        return 72
    elif(character == 'm'):
        return 76
    elif(character == 'q'):
        return 80
    elif(character == 'u'):
        return 84
    elif(character == 'y'):
        return 88
    elif(character == '}'):
        return 92
    else:
        return 93

cdef inline char chr2phImplicit(char character):
    if(character == 33):
        return 0
    elif(character == 34):
        return 1
    elif(character == 35):
        return 2
    elif(character == 36):
        return 3
    elif(character == 37):
        return 4
    elif(character == 38):
        return 5
    elif(character == 39):
        return 6
    elif(character == 40):
        return 7
    elif(character == 41):
        return 8
    elif(character == 42):
        return 9
    elif(character == 43):
        return 10
    elif(character == 44):
        return 11
    elif(character == 45):
        return 12
    elif(character == 46):
        return 13
    elif(character == 47):
        return 14
    elif(character == 48):
        return 15
    elif(character == 49):
        return 16
    elif(character == 50):
        return 17
    elif(character == 51):
        return 18
    elif(character == 52):
        return 19
    elif(character == 53):
        return 20
    elif(character == 54):
        return 21
    elif(character == 55):
        return 22
    elif(character == 56):
        return 23
    elif(character == 57):
        return 24
    elif(character == 58):
        return 25
    elif(character == 59):
        return 26
    elif(character == 60):
        return 27
    elif(character == 61):
        return 28
    elif(character == 62):
        return 29
    elif(character == 63):
        return 30
    elif(character == 64):
        return 31
    elif(character == 65):
        return 32
    elif(character == 66):
        return 33
    elif(character == 67):
        return 34
    elif(character == 68):
        return 35
    elif(character == 69):
        return 36
    elif(character == 70):
        return 37
    elif(character == 71):
        return 38
    elif(character == 72):
        return 39
    elif(character == 73):
        return 40
    elif(character == 74):
        return 41
    elif(character == 75):
        return 42
    elif(character == 76):
        return 43
    elif(character == 77):
        return 44
    elif(character == 78):
        return 45
    elif(character == 79):
        return 46
    elif(character == 80):
        return 47
    elif(character == 81):
        return 48
    elif(character == 82):
        return 49
    elif(character == 83):
        return 50
    elif(character == 84):
        return 51
    elif(character == 85):
        return 52
    elif(character == 86):
        return 53
    elif(character == 87):
        return 54
    elif(character == 88):
        return 55
    elif(character == 89):
        return 56
    elif(character == 90):
        return 57
    elif(character == 91):
        return 58
    elif(character == 92):
        return 59
    elif(character == 93):
        return 60
    elif(character == 94):
        return 61
    elif(character == 95):
        return 62
    elif(character == 96):
        return 63
    elif(character == 97):
        return 64
    elif(character == 98):
        return 65
    elif(character == 99):
        return 66
    elif(character == 100):
        return 67
    elif(character == 101):
        return 68
    elif(character == 102):
        return 69
    elif(character == 103):
        return 70
    elif(character == 104):
        return 71
    elif(character == 105):
        return 72
    elif(character == 106):
        return 73
    elif(character == 107):
        return 74
    elif(character == 108):
        return 75
    elif(character == 109):
        return 76
    elif(character == 110):
        return 77
    elif(character == 111):
        return 78
    elif(character == 112):
        return 79
    elif(character == 113):
        return 80
    elif(character == 114):
        return 81
    elif(character == 115):
        return 82
    elif(character == 116):
        return 83
    elif(character == 117):
        return 84
    elif(character == 118):
        return 85
    elif(character == 119):
        return 86
    elif(character == 120):
        return 87
    elif(character == 121):
        return 88
    elif(character == 122):
        return 89
    elif(character == 123):
        return 90
    elif(character == 124):
        return 91
    elif(character == 125):
        return 92
    else:
        return 93


cdef inline cystr chrACGNTInline(char character):
    if(character == 65):
        return "A"
    elif(character == 67):
        return "C"
    elif(character == 71):
        return "G"
    elif(character == 84):
        return "T"
    elif(character == 78):
        return "N"
    else:
        return "X"


cdef inline cystr chrInline(char character):
    if(character == 33):
            return '!'
    elif(character == 34):
            return '"'
    elif(character == 35):
            return '#'
    elif(character == 36):
            return '$'
    elif(character == 37):
            return '%'
    elif(character == 38):
            return '&'
    elif(character == 39):
            return '\''
    elif(character == 40):
            return '('
    elif(character == 41):
            return ')'
    elif(character == 42):
            return '*'
    elif(character == 43):
            return '+'
    elif(character == 44):
            return ','
    elif(character == 45):
            return '-'
    elif(character == 46):
            return '.'
    elif(character == 47):
            return '/'
    elif(character == 48):
            return '0'
    elif(character == 49):
            return '1'
    elif(character == 50):
            return '2'
    elif(character == 51):
            return '3'
    elif(character == 52):
            return '4'
    elif(character == 53):
            return '5'
    elif(character == 54):
            return '6'
    elif(character == 55):
            return '7'
    elif(character == 56):
            return '8'
    elif(character == 57):
            return '9'
    elif(character == 58):
            return ':'
    elif(character == 59):
            return ';'
    elif(character == 60):
            return '<'
    elif(character == 61):
            return '='
    elif(character == 62):
            return '>'
    elif(character == 63):
            return '?'
    elif(character == 64):
            return '@'
    elif(character == 65):
            return 'A'
    elif(character == 66):
            return 'B'
    elif(character == 67):
            return 'C'
    elif(character == 68):
            return 'D'
    elif(character == 69):
            return 'E'
    elif(character == 70):
            return 'F'
    elif(character == 71):
            return 'G'
    elif(character == 72):
            return 'H'
    elif(character == 73):
            return 'I'
    elif(character == 74):
            return 'J'
    elif(character == 75):
            return 'K'
    elif(character == 76):
            return 'L'
    elif(character == 77):
            return 'M'
    elif(character == 78):
            return 'N'
    elif(character == 79):
            return 'O'
    elif(character == 80):
            return 'P'
    elif(character == 81):
            return 'Q'
    elif(character == 82):
            return 'R'
    elif(character == 83):
            return 'S'
    elif(character == 84):
            return 'T'
    elif(character == 85):
            return 'U'
    elif(character == 86):
            return 'V'
    elif(character == 87):
            return 'W'
    elif(character == 88):
            return 'X'
    elif(character == 89):
            return 'Y'
    elif(character == 90):
            return 'Z'
    elif(character == 91):
            return '['
    elif(character == 92):
            return '\\'
    elif(character == 93):
            return ']'
    elif(character == 94):
            return '^'
    elif(character == 95):
            return '_'
    elif(character == 96):
            return '`'
    elif(character == 97):
            return 'a'
    elif(character == 98):
            return 'b'
    elif(character == 99):
            return 'c'
    elif(character == 100):
            return 'd'
    elif(character == 101):
            return 'e'
    elif(character == 102):
            return 'f'
    elif(character == 103):
            return 'g'
    elif(character == 104):
            return 'h'
    elif(character == 105):
            return 'i'
    elif(character == 106):
            return 'j'
    elif(character == 107):
            return 'k'
    elif(character == 108):
            return 'l'
    elif(character == 109):
            return 'm'
    elif(character == 110):
            return 'n'
    elif(character == 111):
            return 'o'
    elif(character == 112):
            return 'p'
    elif(character == 113):
            return 'q'
    elif(character == 114):
            return 'r'
    elif(character == 115):
            return 's'
    elif(character == 116):
            return 't'
    elif(character == 117):
            return 'u'
    elif(character == 118):
            return 'v'
    elif(character == 119):
            return 'w'
    elif(character == 120):
            return 'x'
    elif(character == 121):
            return 'y'
    elif(character == 122):
            return 'z'
    elif(character == 123):
            return '{'
    elif(character == 124):
            return '|'
    elif(character == 125):
            return '}'
    else:
            return '~'


cdef inline char CigarOpToCigarChar(char character):
    if(character == 0):
        return 77
    elif(character == 1):
        return 73
    elif(character == 2):
        return 68
    elif(character == 3):
        return 83
    elif(character == 4):
        return 78
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
        return 4
    elif(character == 80):
        return 6
    elif(character == 72):
        return 5
    elif(character == 83):
        return 3
    elif(character == 88):
        return 8
    elif(character == 68):
        return 2
    else:
        return 7


cdef inline char ChrToRefIDInline(cystr contig):
    if(contig == 'GL000224.1'):
        return 74
    elif(contig == 'GL000219.1'):
        return 73
    elif(contig == 'GL000205.1'):
        return 72
    elif(contig == 'GL000204.1'):
        return 55
    elif(contig == 'GL000228.1'):
        return 60
    elif(contig == 'GL000198.1'):
        return 56
    elif(contig == 'GL000241.1'):
        return 49
    elif(contig == 'GL000208.1'):
        return 57
    elif(contig == 'GL000212.1'):
        return 77
    elif(contig == 'GL000193.1'):
        return 80
    elif(contig == 'GL000210.1'):
        return 29
    elif(contig == 'GL000239.1'):
        return 30
    elif(contig == 'GL000191.1'):
        return 58
    elif(contig == 'GL000213.1'):
        return 66
    elif(contig == 'GL000220.1'):
        return 65
    elif(contig == 'GL000238.1'):
        return 42
    elif(contig == 'GL000214.1'):
        return 61
    elif(contig == 'GL000194.1'):
        return 81
    elif(contig == 'GL000216.1'):
        return 70
    elif(contig == 'GL000237.1'):
        return 53
    elif(contig == 'GL000243.1'):
        return 50
    elif(contig == 'GL000234.1'):
        return 44
    elif(contig == '20'):
        return 19
    elif(contig == '21'):
        return 20
    elif(contig == '22'):
        return 21
    elif(contig == 'GL000215.1'):
        return 71
    elif(contig == 'GL000196.1'):
        return 39
    elif(contig == 'GL000227.1'):
        return 59
    elif(contig == 'GL000248.1'):
        return 40
    elif(contig == 'GL000197.1'):
        return 35
    elif(contig == 'GL000235.1'):
        return 31
    elif(contig == 'GL000249.1'):
        return 38
    elif(contig == '1'):
        return 0
    elif(contig == '3'):
        return 2
    elif(contig == '2'):
        return 1
    elif(contig == '5'):
        return 4
    elif(contig == '4'):
        return 3
    elif(contig == '7'):
        return 6
    elif(contig == '6'):
        return 5
    elif(contig == '9'):
        return 8
    elif(contig == '8'):
        return 7
    elif(contig == 'GL000199.1'):
        return 68
    elif(contig == 'GL000232.1'):
        return 45
    elif(contig == 'GL000242.1'):
        return 51
    elif(contig == 'GL000236.1'):
        return 48
    elif(contig == 'GL000209.1'):
        return 63
    elif(contig == 'GL000246.1'):
        return 37
    elif(contig == 'GL000244.1'):
        return 41
    elif(contig == 'GL000221.1'):
        return 62
    elif(contig == 'GL000245.1'):
        return 34
    elif(contig == 'GL000203.1'):
        return 36
    elif(contig == 'GL000195.1'):
        return 76
    elif(contig == 'GL000229.1'):
        return 27
    elif(contig == 'GL000226.1'):
        return 26
    elif(contig == 'GL000201.1'):
        return 32
    elif(contig == 'GL000247.1'):
        return 33
    elif(contig == 'Y'):
        return 23
    elif(contig == 'X'):
        return 22
    elif(contig == 'GL000222.1'):
        return 78
    elif(contig == 'GL000202.1'):
        return 43
    elif(contig == 'GL000211.1'):
        return 67
    elif(contig == '11'):
        return 10
    elif(contig == '10'):
        return 9
    elif(contig == '13'):
        return 12
    elif(contig == '12'):
        return 11
    elif(contig == '15'):
        return 14
    elif(contig == '14'):
        return 13
    elif(contig == '17'):
        return 16
    elif(contig == '16'):
        return 15
    elif(contig == '19'):
        return 18
    elif(contig == '18'):
        return 17
    elif(contig == 'GL000231.1'):
        return 28
    elif(contig == 'GL000233.1'):
        return 54
    elif(contig == 'GL000240.1'):
        return 47
    elif(contig == 'GL000200.1'):
        return 79
    elif(contig == 'GL000230.1'):
        return 52
    elif(contig == 'GL000217.1'):
        return 69
    elif(contig == 'MT'):
        return 24
    elif(contig == 'GL000223.1'):
        return 75
    elif(contig == 'GL000225.1'):
        return 82
    elif(contig == 'GL000206.1'):
        return 46
    elif(contig == 'GL000218.1'):
        return 64
    elif(contig == 'GL000207.1'):
        return 25
    else:
        return 83
