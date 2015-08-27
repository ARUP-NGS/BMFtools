#include "dmp.h"

//Function definitions
double_t igamc_pvalues(int num_pvalues, double_t x);
double_t INV_CHI2_FROM_LOG10(int32_t log10int);

int main(int argc, char* argv[]) {
    int n;
    double f;
    n = atoi(argv[1]);
    f = atof(argv[2]);
    fprintf(stderr, "IGAMC_PVALUES: %f", igamc_pvalues((double_t)n, (double_t)f));
    return 0;
}
