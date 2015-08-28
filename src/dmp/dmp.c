#include <getopt.h>
#include <quadmath.h>
#include "dmp.h"

//Function definitions
float128_t igamc_pvalues(int num_pvalues, float128_t x);
KingFisher_t init_kf(int max, size_t readlen);
void pushback_kseq(KingFisher_t *fisher, kseq_t *seq);

void print_usage() {
	fprintf(stderr, "This isn't ready to do anything yet. Oops.\n");
}

void print_opt_err(char *argv[], char optarg[]) {
	fprintf(stderr, "Invalid argument %s. See usage.\n", optarg);
	print_usage();
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
#if !NDEBUG
    fprintf(stderr, "IGAMC_PVALUES: %.30Lf\n", (long double)igamc_pvalues((float128_t)n, (float128_t)f));
    fprintf(stderr, "degrees of freedom:  %i. Chi2 sum: %f. Chi2 sum as a string: %s.\n", n, f, argv[2]);
    float128_t *arr = (float128_t *)calloc(1, sizeof(float128_t));
    int wtf = arr[0];
    fprintf(stderr, "A zeroed 128-bit float has value %i\n", wtf);
#endif
    return 0;
}
