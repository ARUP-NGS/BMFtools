#ifndef BAM_PR_H
#define BAM_PR_H
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <zlib.h>
#include <tgmath.h>
#include "htslib/sam.h"
#include "include/sam_opts.h"
#include "include/bam.h" // for bam_get_library
#include "include/igamc_cephes.h" /// for igamc
#include "dlib/cstr_util.h"
#include "dlib/sort_util.h"
#include "dlib/bam_util.h"

#define STACK_START 128
#define SEQBUF_SIZE 300

#define seq2buf(buf, seq, len) \
	do {\
		uint64_t i_##seq;\
		for(i_##seq = 0; i_##seq < (len >> 1); ++i_##seq) {\
			buf[i_##seq] = seq_nt16_str[bam_seqi(seq, i_##seq)];\
			buf[len - i_##seq - 1] = seq_nt16_str[bam_seqi(seq, len - i_##seq - 1)];\
		}\
		if(len&1) buf[i_##seq] = seq_nt16_str[bam_seqi(seq, i_##seq)];\
		buf[len] = '\0';\
	} while(0)



typedef bam1_t *bam1_p;

typedef struct {
	int n, max;
	bam1_t **a;
} tmp_stack_t;

const int8_t seq_comp_table[16] = { 0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15 };

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

typedef int (*stack_fn)(bam1_t *b, bam1_t *p);

CONST static inline int same_stack_pos_se(bam1_t *b, bam1_t *p)
{
	return bmfsort_se_key(b) == bmfsort_se_key(p);
}

CONST static inline int same_stack_ucs_se(bam1_t *b, bam1_t *p)
{
	return ucs_se_sort_key(b) == ucs_se_sort_key(p);
}

CONST static inline int same_stack_pos(bam1_t *b, bam1_t *p)
{
	return (bmfsort_core_key(b) == bmfsort_core_key(p) &&
			bmfsort_mate_key(b) == bmfsort_mate_key(p));
}

CONST static inline int same_stack_ucs(bam1_t *b, bam1_t *p)
{
#if !NDEBUG
	if(!p) {
		fprintf(stderr, "Later bam record null. Abort!\n");
		exit(EXIT_FAILURE);
	}
	if(!b) {
		fprintf(stderr, "First bam record null. Abort!\n");
		exit(EXIT_FAILURE);
	}
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
    stack_fn fn;
} pr_settings_t;


CONST static inline int read_pass_hd(bam1_t *b, bam1_t *p, const int lim)
{
	const uint8_t *const bseq = bam_get_seq(b);
	const uint8_t *const pseq = bam_get_seq(p);
	int hd = 0;
	for(int i = 0; i < b->core.l_qseq; ++i) {
		const uint8_t bc = bam_seqi(bseq, i);
		const uint8_t pc = bam_seqi(pseq, i);
		if(bc != pc && bc != HTS_N && pc != HTS_N)
			if(++hd > lim)
				return 0;
	}
}

int bam_rsq(int argc, char *argv[]);
void bam2ffq(bam1_t *b, FILE *fp);
void write_stack(tmp_stack_t *stack, pr_settings_t *settings);



#endif /* BAM_PR_H */
