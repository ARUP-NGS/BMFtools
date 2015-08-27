#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

#include "dmp_utils.h"


#define LOG10E_X5_INV 0.4605170185988091368035982909368728415202202977257545952066655801935145219354704960471994410179196596683935568084572497266819050930165613513332L
//Multiply a phred score by this to convert a -10log_10(x) to a -2log_e(x)
//such as in the following macro:
#define LOG10_TO_CHI2(x) (x) * LOG10_X5_INV

#define INV_CHI2_FROM_LOG10(log10int) -2 * log(1 - pow(10, log10int))
/*
 * Equivalent to the following, but type-general:
inline longdouble_t INV_CHI2_FROM_LOG10(int32_t log10int)
{
    return -2 * log(1 - pow(10, log10int));
}
*/

extern long double igamcl(long double a, long double x);

// Converts a 
inline longdouble_t igamc_pvalues(int num_pvalues, longdouble_t x)
{
    if(x < 0) {
        return 1.0;
    }
    else {
#if !NDEBUG
    	fprintf(stderr, "Now calling igamcl.\n");
#endif
        return igamcl(num_pvalues * 1., x / 2.0);
    }
}
