#ifndef BAM_PR_H
#define BAM_PR_H
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <zlib.h>
#include <unistd.h>
#include <tgmath.h>
#include <inttypes.h>
#include "htslib/sam.h"
#include "sam_opts.h"
#include "bam.h" // for bam_get_library
#include "bam_rescue.h"
#include "igamc_cephes.c"
#include "cstr_util.h"

#define STACK_START (1 << 5)

//Multiply a phred score by this to convert a -10log_10(x) to a -2log_e(x)
#define LOG10E_X5_INV 0.460517018598809136803598290936872841520220297725754595206665580193514521935470496
#define LOG10E_X5_1_2 0.230258509299404568401799145468436420760110148862877297603332790096757260967735248
//such as in the following macro:
#define LOG10_TO_CHI2(x) (x) * LOG10E_X5_INV


#ifndef bam_is_r1
#define bam_is_r1(b) (!!((b)->core.flag&BAM_FREAD1))
#endif

#ifndef INC_TAG
#define INC_TAG(p, b, key) *(int *)(bam_aux_get(p, key) + 1) += *(int *)(bam_aux_get(b, key) + 1);
#endif

#ifndef bam_is_r2
#define bam_is_r2(b) (!!((b)->core.flag&BAM_FREAD2))
#endif


#ifndef ucs_sort_mate_key
#define ucs_sort_mate_key(a) (uint64_t)((uint64_t)a->core.mtid<<32|(bam_aux2i(bam_aux_get(b, "MU")) + 1)<<2|(!!(a->core.flag & BAM_FMREVERSE)))
#endif

#ifndef ucs_sort_core_key
#define ucs_sort_core_key(a) (uint64_t)((uint64_t)a->core.tid<<32|(bam_aux2i(bam_aux_get(b, "SU"))+1)<<2|bam_is_rev(a)<<1|bam_is_r1(a))
#endif

#ifndef bam_sort_core_key
#define bam_sort_core_key(a) (uint64_t)((uint64_t)a->core.tid<<32|(a->core.pos+1)<<2|bam_is_rev(a)<<1|bam_is_r1(a))
#endif

#ifndef bam_sort_mate_key
#define bam_sort_mate_key(a) (uint64_t)((uint64_t)a->core.mtid<<32|(a->core.mpos+1)<<1|(!!(a->core.flag & BAM_FMREVERSE)))
#endif

#ifndef AVG_LOG_TO_CHI2
#define AVG_LOG_TO_CHI2(x, y) (x + y) * LOG10E_X5_1_2
#endif

#ifndef HTS_A
#define HTS_A 1
#endif

#ifndef HTS_C
#define HTS_C 2
#endif

#ifndef HTS_G
#define HTS_G 4
#endif

#ifndef HTS_T
#define HTS_T 8
#endif

#ifndef HTS_N
#define HTS_N 15
#endif

#define SEQBUF_SIZE 300

#define seq2buf(buf, seq, len) \
	uint64_t i_;\
	for(i_ = 0; i_ < len; ++i_) {\
		buf[i_] = seq_nt16_str[bam_seqi(seq, i_)];\
		buf[len - i_ - 1] = seq_nt16_str[bam_seqi(seq, len - i_ - 1)];\
	}\
	if(len&1)\
		buf[i_] = seq_nt16_str[bam_seqi(seq, i_)];\
	buf[len] = '\0'


//#define set_base(pSeq, bSeq, i) (pSeq)[(i)>>1] = (i % 2) ? (bam_seqi((bSeq), i) | ((pSeq)[(i)>>1] & 0xf0U)): ((bam_seqi((bSeq), i) << 4) | ((pSeq)[(i)>>1] & 0xfU))
#define set_base(pSeq, bSeq, i) (pSeq)[(i)>>1] = ((bam_seqi((bSeq), i) << ((!(i % 2)) << 2)) | (((pSeq)[(i)>>1]) & (0xfU << ((!(i % 2)) << 2))))


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

static inline void resize_stack(tmp_stack_t *stack, size_t n) {
	if(!stack->a) {
		stack->a = (bam1_t **)malloc(sizeof(bam1_t *) * n);
		memset(stack->a, 0, sizeof(bam1_t *) * n);
		stack->max = n;
		return;
	}
	if(n > stack->max) {
		stack->max = n;
		stack->a = (bam1_t **)realloc(stack->a, sizeof(bam1_t *) * n);
#if !NDEBUG
		if(!stack->a) {
			fprintf(stderr, "Failed to reallocate memory for %i bam1_t * objects. Abort!\n", stack->max);
			exit(EXIT_FAILURE);
		}
#endif
	}
	else if(n < stack->n){
		for(uint64_t i = stack->n; i > n; --i) {
			bam_destroy1(stack->a[i]);
		}
		stack->max = n;
		stack->a = (bam1_t **)realloc(stack->a, sizeof(bam1_t *) * n);
	}
}

enum cmpkey {
	POS,
	UCS
};

static inline int same_stack_pos(bam1_t *b, bam1_t *p) {
	return (bam_sort_core_key(b) == bam_sort_core_key(p) &&
			bam_sort_mate_key(b) == bam_sort_mate_key(p));
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
	int annealed; // Set to true to check a reversed barcode for a mismatch limit.
	int is_se; // Is single-end
	bam_hdr_t *hdr; // BAM header

} pr_settings_t;

static inline int pvalue_to_phred(double pvalue)
{
	return (int)(-10 * log10(pvalue) + 0.5); // Add 0.5 to round up
}

static inline int disc_pvalues(double pv_better, double pv_worse)
{
	return pvalue_to_phred(igamc(2., LOG10_TO_CHI2(pv_better - (10. * log10(1 - pow(10., (pv_worse * 0.1)))))));
}


// Converts a chi2 sum into a p value.
static inline double igamc_pvalues(int num_pvalues, double x)
{
	if(x < 0) {
		return 1.0;
	}
	else {
		return igamc((double)num_pvalues, x / 2.0);
	}
}

static inline void *array_tag(bam1_t *b, const char *tag) {
	uint8_t *data = bam_aux_get(b, tag);
#if !NDEBUG
	if(*data++ != 'B') {
		fprintf(stderr, "This is not an array tag. Abort mission! (%c)\n", *data);
		exit(EXIT_FAILURE);
	}
	char typecode = *data++;
	int n = *((int *)data);
	return data ? (void *)(data + sizeof(int)): NULL; // Offset by 1 to skip typecode, 2 to skip array length.
#else
	return data ? (void *)(data + sizeof(int)): NULL;
#endif
}


static inline void bam2ffq(bam1_t *b, FILE *fp)
{
	char comment[3000] = "";
	uint32_t *pv = (uint32_t *)array_tag(b, (char *)"PV");
	uint32_t *fa = (uint32_t *)array_tag(b, (char *)"FA");
	append_csv_buffer(b->core.l_qseq, (int *)pv, comment, (char *)"PV:I:B");
	strcat(comment, "\t");
	append_csv_buffer(b->core.l_qseq, (int *)fa, comment, (char *)"FA:I:B");
	append_int_tag(comment, (char *)"FM", bam_aux2i(bam_aux_get(b, (char *)"FM")));
	append_int_tag(comment, (char *)"RV", bam_aux2i(bam_aux_get(b, (char *)"RV")));
	append_int_tag(comment, (char *)"FP", bam_aux2i(bam_aux_get(b, (char *)"FP")));
	if(!(b->core.flag & BAM_FREVERSE)) {
		fprintf(fp, "Writing the read with current direction.\n");
		uint8_t *bseq = bam_get_seq(b);
		uint8_t *bqual = bam_get_qual(b);
		char seqbuf[SEQBUF_SIZE];
		seq2buf(seqbuf, bseq, b->core.l_qseq);
		fprintf(fp, "@%s %s\n%s\n+\n", bam_get_qname(b), comment, seqbuf);
		for(int i = 0; i < b->core.l_qseq; ++i)
			seqbuf[i] = 33 + bqual[i];
		fprintf(fp, "%s\n", seqbuf);
		return;
	}
	fprintf(stderr, "Writing the read by reverse-complementing.\n");
	uint8_t *bseq = bam_get_seq(b);
	uint8_t *bqual = bam_get_qual(b);
	char seqbuf[SEQBUF_SIZE];
	int i;
	for(i = 0; i < (b->core.l_qseq >> 1); ++i) {
		const int8_t tmp = bam_seqi(bseq, i);
		seqbuf[i] = nuc_cmpl(seq_nt16_str[bam_seqi(bseq, b->core.l_qseq - i - 1)]);
		seqbuf[b->core.l_qseq - i - 1] = nuc_cmpl(seq_nt16_str[tmp]);
	}
	if(b->core.l_qseq & 1)
		seqbuf[i] = nuc_cmpl(seq_nt16_str[bam_seqi(bseq, i)]);
#if !NDEBUG
	fprintf(stderr, "Now the sequence is %s.\n", seqbuf);
#endif
	fprintf(fp, "Writing the read by reverse-complementing.\n");
	fprintf(fp, "@%s %s\n%s\n+\n", (char *)bam_get_qname(b), comment, seqbuf);
	for(uint64_t j = 0; j < b->core.l_qseq; ++j)
		seqbuf[j] = bqual[b->core.l_qseq - j - 1] + 33;
	fprintf(fp, "%s\n", seqbuf);
	return;
}


static inline void write_stack(tmp_stack_t *stack, pr_settings_t *settings)
{
	for(int i = 0; i < stack->n; ++i) {
		if(stack->a[i]) {
			if(bam_aux_get(stack->a[i], "RA")) { /* If RA tag is present, IE, if it was merged.*/
				//fprintf(stderr, "This should be writing the record to the fastq handle.\n");
				bam2ffq(stack->a[i], settings->fqh);
			}
			else {
				//fprintf(stderr, "This should be writing the record to the bam handle.\n");
				sam_write1(settings->out, settings->hdr, stack->a[i]);
			}
			bam_destroy1(stack->a[i]);
			stack->a[i] = NULL;
		}
	}
}

int bam_pr(int argc, char *argv[]);



#endif /* BAM_PR_H */
