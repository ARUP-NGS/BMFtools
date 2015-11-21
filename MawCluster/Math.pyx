# cython: boundscheck=False, wraparound=False, initializedcheck=False
"""
Port of gamma functions for calculating Fisher scores.
"""

DEF LOG10E_X5_INV = 0.46051701859880917
# Multiply a phred score by this to convert

cdef inline double_t CHI2_FROM_PHRED(int32_t phredInt) nogil:
    return phredInt * LOG10E_X5_INV


cdef inline double_t INV_CHI2_FROM_PHRED(int32_t phredInt) nogil:
    return -2 * c_log(1 - c_pow(10, (phredInt * -1.)))


cdef int8_t * arrmax(int32_t * quals, int8_t * ret,
                     size_t nRecs, size_t rLen) nogil:
    """
    :param quals: [int32_t **/arg] two-dimensional array of
    32-bit quality scores. First dimension, read, second dimension, cycle
    :param ret: [int32_t */arg] one-dimensional array
    into which to copy the results. This should be zero'd.
    :param nRecs: [size_t/arg] number of records
    :param rLen: [size_t/arg] length of the reads.
    :return: void
    """
    cdef size_t index1, index2
    for index1 in range(nRecs):
        for index2 in range(rLen):
            if quals[index2 + index1 * rLen] > ret[index2]:
                ret[index2] = quals[index2 + index1 * rLen] + 33
    return ret
