#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "cephes/igam.c"

#define LOG10E_X5_INV = 0.46051701859880917
//Multiply a phred score by this to convert a -10log_10(x) to a -2log_e(x)
//such as in the following macro:
#define LOG10_TO_CHI2(x) (x) * LOG10_X5_INV

typedef double double_t;

// Converts a 
inline double_t igamc_pvalues(int num_pvalues, double_t x)
{
    if(x < 0) {
        return 1.0;
    }
    else {
        return igamc(num_pvalues * 1., x / 2.0);
    }
}

inline double_t INV_CHI2_FROM_LOG10(int32_t log10int)
{
    return -2 * log(1 - pow(10, log10int));
}
