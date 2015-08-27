/*
 * Utilities for performing molecular demultiplexing
 */

#include "include/kseq.h"
#include <zlib.h>
// Force declaration of all of kseq's types.
KSEQ_INIT(gzFile, gzread)

#define KF_MAX_INCREMENT 20

typedef long double longdouble_t;

/*
 * TODO: KingFisher finishing work.
 * A destructor for KingFisher.
 * Rewrite the cFastFisherFlattening array work in C rather than Cython from MawCluster/BCFastq.pyx.
 * Use that array work to fill in the update_kf method.
 * Ultimately, I would like to only keep the seq strings in memory.
 * Or, alternatively, I could be keeping only counts of nucleotides! Much better way.
 */

typedef struct KingFisher {
	char *barcode;
	char **seqs;
	longdouble_t **chi2sums;
	int length;
	int max;
} KingFisher_t;

inline KingFisher_t *init_kf(int max, size_t readlen) {
	KingFisher_t fisher = {
		.barcode = NULL,
		.seqs = (char *)malloc(max * sizeof(char *)),
		.chi2sums = (longdouble_t **)malloc(sizeof(longdouble_t *) * 4),
		.length = 0,
		.max = max
	}
	return &fisher;
}

inline void update_kf(KingFisher_t *fisher, kseq_t *seq){
	fprintf(stderr, "update_kf for updating KingFisher_t is unimplemented. Abort!\n");
	exit(1);
	return;
}

inline void pushback_kseq(KingFisher_t *fisher, kseq_t *seq) {
#if !NDEBUG
	if(fisher->length > fisher->max) {
		fprintf(stderr, "KingFisher length %i is greater than its max %i. ???\n", fisher->length, fisher->max);
		exit(1);
	}
#endif

	if(fisher->max == fisher->length) {
		fisher->max += KF_MAX_INCREMENT;
		fisher->seqs = (char **)realloc(fisher->seqs, fisher->max * sizeof(char *));
	}
	fisher->seqs[length] = strdup(seq->seq.s);
	fisher->length++;
	update_kf(fisher, seq);
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
