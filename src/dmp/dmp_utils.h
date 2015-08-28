/*
 * Utilities for performing molecular demultiplexing
 */

#include <zlib.h>
#include <stddef.h>
#include <quadmath.h>
#include "include/kseq.h"
// Force declaration of all of kseq's types.
KSEQ_INIT(gzFile, gzread)

typedef __float128 float128_t;

/*
 * TODO: KingFisher finishing work.
 * A destructor for KingFisher.
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
		nuc_counts[i] = (int *)calloc(5, sizeof(int));
		chi2sums[i] = (float128_t *)calloc(4, sizeof(float128_t));
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
/*
inline void update_nuc_counts(KingFisher_t *fisher, kseq_t *seq){

	fprintf(stderr, "update_kf for updating KingFisher_t is unimplemented. Abort!\n");
	exit(1);
	return;
}
*/

inline void pushback_kseq(KingFisher_t *fisher, kseq_t *seq) {
	for(int i = 0; i < fisher->readlen; i++) {

	}
	fisher->length++; // Increment
	return;
}

/*
 *
 * I'm not sure if I'll be using this directly.
typedef struct SeqQual {
	char *seq;
	char *qual;
	int length;
} SeqQual_t;

#define INIT_SEQ_QUAL(var, kseq_ptr) SeqQual_t var = {           \
		.seq = strdup(kseq_ptr->s),                              \
		.qual = strdup(kseq_ptr->s),                             \
		.length = kseq_ptr->l                                    \
		};


#define DESTROY_SEQ_QUAL(var) free(var.seq);                     \
		free(var.qual);                                          \
*/
