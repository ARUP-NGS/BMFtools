#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <stddef.h>
#include <quadmath.h>
#include "include/kseq.h"

// Force declaration of all of kseq's types.
KSEQ_INIT(gzFile, gzread)

typedef __float128 float128_t;


#define LOG10E_X5_INV 0.4605170185988091368035982909368728415202202977257545952066655801935145219354704960471994410179196596683935568084572497266819050930165613513332L
//Multiply a phred score by this to convert a -10log_10(x) to a -2log_e(x)
//such as in the following macro:
#define LOG10_TO_CHI2(x) (x) * LOG10E_X5_INV

#define INV_CHI2_FROM_LOG10(log10int) -2 * log(1 - pow(10, log10int))
/*
 * Equivalent to the following, but type-general:
inline float128_t INV_CHI2_FROM_LOG10(int32_t log10int)
{
    return -2 * log(1 - pow(10, log10int));
}
*/


extern float128_t igamcl(float128_t a, float128_t x);

// Converts a 
inline float128_t igamc_pvalues(int num_pvalues, float128_t x)
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

typedef struct intpair {
    int i1;
    int i2;
} intpair_t;


inline intpair_t NUC_TO_POS(char character) {
    intpair_t ret;
    switch(character) {
        case 'A': ret.i1 = 0; ret.i2 = 0; break;
        case 'C': ret.i1 = 1; ret.i2 = 1; break;
        case 'G': ret.i1 = 2; ret.i2 = 2; break;
        case 'T': ret.i1 = 3; ret.i2 = 3; break;
        default: ret.i1 = 0; ret.i2 = 4; break;
    }
    return ret;
}

/*
 * TODO: KingFisher finishing work.
 * A destructor for KingFisher.
 * FFF core.
 * Output result:
 * Rewrite the cFastFisherFlattening array work in C rather than Cython from MawCluster/BCFastq.pyx.
 * Use that array work to fill in the update_kf method.
 */

typedef struct KingFisher {
    char *barcode; // Barcode for the family
    int **nuc_counts; // Count of nucleotides of this form
    float128_t **chi2sums; // Sums of -2ln(p-value)
    int length; // Number of reads in family
    int readlen; // Length of reads
} KingFisher_t;

inline KingFisher_t init_kf(int max, size_t readlen) {
    int **nuc_counts = (int **)malloc(readlen * sizeof(int *));
    float128_t **chi2sums = (float128_t **)malloc(sizeof(float128_t *) * 4);
    for(int i = 0; i < readlen; i++) {
        nuc_counts[i] = (int *)calloc(5, sizeof(int)); // One each for A, C, G, T, and N
        chi2sums[i] = (float128_t *)calloc(4, sizeof(float128_t)); // One for each nucleotide
    }
    KingFisher_t fisher = {
        .barcode = NULL,
        .nuc_counts = nuc_counts,
        .chi2sums = chi2sums,
        .length = 0,
        .readlen = readlen
    };
    return fisher;
}

inline void destroy_kf(KingFisher_t *kfp) {
    for(int i = 0; i < kfp->readlen; i++) {
        free(kfp->nuc_counts[i]);
        free(kfp->chi2sums[i]);
    }
    free(kfp->barcode);
    free(kfp->nuc_counts);
    free(kfp->chi2sums);
}

inline void clear_kf(KingFisher_t *kfp) {
    for(int i = 0; i < kfp->readlen; i++) {
        memset(kfp->chi2sums[i], 0, 4 * sizeof(float128_t)); // Sets these to 0.
        memset(kfp->nuc_counts[i], 0, 5 * sizeof(int)); // And these.
    }
    kfp->length = 0;
    return;
}

/*
inline void update_nuc_counts(KingFisher_t *fisher, kseq_t *seq){

    fprintf(stderr, "update_kf for updating KingFisher_t is unimplemented. Abort!\n");
    exit(1);
    return;
}
*/

inline void pushback_kseq(KingFisher_t *kfp, kseq_t *seq) {
    intpair_t nuc_indices;
    for(int i = 0; i < kfp->readlen; i++) {
        nuc_indices = NUC_TO_POS((seq->seq.s[i]));
        kfp->nuc_counts[i][nuc_indices.i2] += 1;
        kfp->chi2sums[i][nuc_indices.i1] += LOG10_TO_CHI2((seq->qual.s[i] - 33));
    }
    kfp->length++; // Increment
    return;
}
