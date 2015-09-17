#include <stdint.h>
#include <stdlib.h>
#include <math.h>
typedef double float64_t;
#define square(x) ((x)*(x))

float64_t Hellinger_in_c(float64_t* arr1, float64_t* arr2, size_t length){
    float64_t cumSum = 0.;
    float64_t tmpFloat;
    for(size_t index = 0; index < length; index++){
        tmpFloat = abs(sqrt(arr1[index]) - sqrt(arr2[index]));
        cumSum += square(tmpFloat);
    }
    return sqrt(cumSum) * M_SQRT1_2;
}
/*
 * cdef float64_t cHellingerDistance(float64_t* array1,
                                  float64_t* array2,
                                  size_t length) nogil:
    """
    Calculates the Helling Distance between two discrete probability
    distributions, described (each) as a 1-dimensional array.
    """
    cdef float64_t cumSum = 0.
    cdef float64_t tmpFloat
    cdef size_t tmpInt
    for tmpInt in range(length):
        tmpFloat = c_abs(c_sqrt(array1[tmpInt]) -
                         c_sqrt(array2[tmpInt]))
        cumSum += c_square(tmpFloat)
    return c_sqrt(cumSum) * M_SQRT1_2
 */
