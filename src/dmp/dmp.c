#include "dmp.h"

//Function definitions
longdouble_t igamc_pvalues(int num_pvalues, double_t x);
longdouble_t INV_CHI2_FROM_LOG10(int32_t log10int);

int main(int argc, char* argv[]) {
    int n;
    double f;
    n = atoi(argv[1]);
    f = atof(argv[2]);
    fprintf(stderr, "IGAMC_PVALUES: %f", igamc_pvalues((longdouble_t)n, (longdouble_t)f));
    return 0;
}
