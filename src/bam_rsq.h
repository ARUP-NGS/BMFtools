#ifndef BAM_PR_H
#define BAM_PR_H
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <zlib.h>
#include <unistd.h>
#include <tgmath.h>
#include <inttypes.h>
#include <assert.h>
#include "htslib/sam.h"
#include "sam_opts.h"
#include "bam.h" // for bam_get_library
#include "igamc_cephes.h"
#include "cstr_util.h"
#include "sort_util.h"
#include "pair_util.h"
#include "compiler_util.h"
//#include "bmfsort.h"

#define STACK_START (1 << 5)



#ifndef INC_TAG
#define INC_TAG(p, b, key) *(int *)(bam_aux_get(p, key) + 1) += *(int *)(bam_aux_get(b, key) + 1);
#endif


#define SEQBUF_SIZE 300

#define seq2buf(buf, seq, len) \
	uint64_t i_##seq;\
	for(i_##seq = 0; i_##seq < (len >> 1); ++i_##seq) {\
		buf[i_##seq] = seq_nt16_str[bam_seqi(seq, i_##seq)];\
		buf[len - i_##seq - 1] = seq_nt16_str[bam_seqi(seq, len - i_##seq - 1)];\
	}\
	if(len&1)\
		buf[i_##seq] = seq_nt16_str[bam_seqi(seq, i_##seq)];\
	buf[len] = '\0'

/*
#define set_base(pSeq, bSeq, i) (pSeq)[(i)>>1] = (i % 2) ? \
		(bam_seqi((bSeq), i) | ((pSeq)[(i)>>1] & 0xf0U)): \
		((bam_seqi((bSeq), i) << 4) | ((pSeq)[(i)>>1] & 0xfU))
*/
#define set_base(pSeq, bSeq, i) (pSeq)[(i)>>1] = ((bam_seqi(bSeq, i) << (((~i) & 1) << 2)) | (((pSeq)[(i)>>1]) & (0xf0U >> (((~i) & 1) << 2))))
#define n_base(pSeq, i) pSeq[(i)>>1] |= (0xf << ((~(i) & 1) << 2));

#define check_fa(arr, fm, len) \
		do {\
		for(int i##arr = 0; i##arr < len; ++i##arr) {\
			if(arr[i##arr] > fm){\
				fprintf(stderr, "FAIL. %" PRIu32 " arr value greater than FM %" PRIu32 ".\n", arr[i##arr], fm);\
				exit(EXIT_FAILURE);\
			}\
		}\
		} while(0)


typedef bam1_t *bam1_p;

typedef struct {
	int n, max;
	bam1_t **a;
} tmp_stack_t;

const int8_t seq_comp_table[16] = { 0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15 };

extern double igamc(double a, double x);
extern inline uint32_t pvalue_to_phred(double pvalue);
extern inline double igamc_pvalues(int num_pvalues, double x);

static inline void stack_insert(tmp_stack_t *stack, bam1_t *b)
{
	if (stack->n == stack->max) {
	    stack->max = stack->max? stack->max<<1 : 0x10000;
	    stack->a = (bam1_t**)realloc(stack->a, sizeof(bam1_t*) * stack->max);
	}
	stack->a[stack->n++] = bam_dup1(b);
}

void resize_stack(tmp_stack_t *stack, size_t n);


enum cmpkey {
	POSITION,
	UNCLIPPED
};

static inline int same_stack_pos(bam1_t *b, bam1_t *p) {
	return (bmfsort_core_key(b) == bmfsort_core_key(p) &&
			bmfsort_mate_key(b) == bmfsort_mate_key(p));
}

static inline int same_stack_ucs(bam1_t *b, bam1_t *p) {
#if !NDEBUG
	if(!p) {
		fprintf(stderr, "Later bam record not set. Abort!\n");
		exit(EXIT_FAILURE);
	}
	if(!b) {
		fprintf(stderr, "First bam record not set. Abort!\n");
		exit(EXIT_FAILURE);
	}
	//fprintf(stderr, "Core key 1: %" PRIu64 ". Core key 2: %" PRIu64 ".\n", ucs_sort_core_key(p), ucs_sort_core_key(b));
	//fprintf(stderr, "Mate key 1: %" PRIu64 ". Mate key 2: %" PRIu64 ".\n", ucs_sort_mate_key(p), ucs_sort_mate_key(b));
#endif
	return (ucs_sort_core_key(b) == ucs_sort_core_key(p) &&
			ucs_sort_mate_key(b) == ucs_sort_mate_key(p));
}


typedef struct pr_settings {
	FILE *fqh;
	samFile *in;
	samFile *out;
	int cmpkey; // 0 for pos, 1 for unclipped start position
	int mmlim; // Mismatch failure threshold.
	int realign_unchanged; // Set to true to realign unchanged reads.
	int is_se; // Is single-end
	bam_hdr_t *hdr; // BAM header
} pr_settings_t;


static inline void *array_tag(bam1_t *b, const char *tag) {
	uint8_t *data = bam_aux_get(b, tag);
#if !NDEBUG
	if(*data++ != 'B') {
		fprintf(stderr, "This is not an array tag. Abort mission! (%c)\n", *data);
		exit(EXIT_FAILURE);
	}
	char typecode = *data++;
	//int n = *((int *)data);
	return data ? (void *)(data + sizeof(int)): NULL; // Offset by 1 to skip typecode, 2 to skip array length.
#else
	return data ? (void *)(data + sizeof(int) + 2): NULL;
#endif
}

int bam_rsq(int argc, char *argv[]);
void bam2ffq(bam1_t *b, FILE *fp);
void write_stack(tmp_stack_t *stack, pr_settings_t *settings);



#endif /* BAM_PR_H */
