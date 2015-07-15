# Math functions

ctypedef double double_t
from libc.math cimport exp, INFINITY, M_PI as PI, floor, sin, log, fabs
from libc.stdlib cimport malloc, free, realloc


# CONSTANTS TAKEN FROM CEPHES under UNK setting.
cdef int sgngam = 0
cdef char sizedouble = 8
cdef double_t LOGPI = 1.14472988584940017414
cdef double_t MAXLGM = 2.035093e36
cdef double_t LS2PI  =  0.91893853320467274178
cdef double_t MAXLOG = 8.8029691931113054295988E1
cdef double_t MACHEP =  1.11022302462515654042E-16
cdef double_t big = 4.503599627370496e15
cdef double_t biginv =  2.22044604925031308085e-16


cdef double_t[5] A = [8.11614167470508450300E-4,
                      -5.95061904284301438324E-4,
                      7.93650340457716943945E-4,
                      -2.77777777730099687205E-3,
                      8.33333333333331927722E-2]
cdef double_t[6] B = [-1.37825152569120859100E3,
                     -3.88016315134637840924E4,
                     -3.31612992738871184744E5,
                     -1.16237097492762307383E6,
                     -1.72173700820839662146E6,
                     -8.53555664245765465627E5]
cdef double_t[6] C = [-3.51815701436523470549E2, -1.70642106651881159223E4,
                      -2.20528590553854454839E5, -1.13933444367982507207E6,
                      -2.53252307177582951285E6, -2.01889141433532773231E6]

# FUNCTIONS
cdef inline double_t igamc(double_t a, double_t x) nogil