#include "dmp.h"
#include <getopt.h>

//Function definitions
longdouble_t igamc_pvalues(int num_pvalues, longdouble_t x);

void print_opt_err(char *argv[], char optarg[]) {
	fprintf(stderr, "Invalid argument %s. See usage.\n", optarg);
	exit(1);
}

int main(int argc, char* argv[]) {
    int n = 0;
    double f = 1337.;
    int c;
    while ((c = getopt(argc, argv, "n:f:")) > -1) {
        switch(c) {
            case 'n': n = atoi(optarg); break;
            case 'f': f = atof(optarg); break;
            default: print_opt_err(argv, optarg);
        }
    }
    fprintf(stderr, "IGAMC_PVALUES: %.12Lf\n", igamc_pvalues((longdouble_t)n, (longdouble_t)f));
    fprintf(stderr, "degrees of freedom:  %i. Chi2 sum: %f. Chi2 sum as a string: %s.\n", n, f, argv[2]);
    return 0;
}
