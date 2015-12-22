#ifndef MSEQ_H
#define MSEQ_H
#include <stdio.h>
#include <stdint.h>
#include <zlib.h>
#include "cstr_util.h"
#include "compiler_util.h"
#include "mem_util.h"
#include "rescaler.h"


#define MAX_BARCODE_LENGTH 30 // Maximum expected inline barcode

// Struct definitions


typedef struct mseq {
	char name[100];
	char comment[2000];
	char seq[200];
	char qual[200];
	char barcode[MAX_BARCODE_LENGTH + 1];
	int l;
	int blen;
} mseq_t;

typedef struct tmp_mseq {
	char *tmp_seq;
	char *tmp_qual;
	char *tmp_barcode;
	int readlen;
	int blen;
} tmp_mseq_t;

// KSEQ Utilities
#include "kseq.h"

#ifndef KSEQ_DEC_GZ
#define KSEQ_DEC_GZ
KSEQ_INIT(gzFile, gzread)
#endif

// Functions

static inline int set_barcode(kseq_t *seq1, kseq_t *seq2, char *barcode, int offset, int blen1_2);
CONST static inline char rescale_qscore(int readnum, char qscore, int cycle, char base, int readlen, char *rescaler);
char *mem_view(char *);

// calls mem_view on the comment field of the kseq_t struct.
#define barcode_mem_view(seq) mem_view(seq->comment.s)

/*
 * @func switch_test
 * Decides whether or not to switch read 1 and read 2 based on the barcode sequence.
 * :param: seq1 [kseq_t *] Read 1 kseq object.
 * :param: seq2 [kseq_t *] Read 2 kseq object.
 * :param:
 */
CONST static inline int switch_test(kseq_t *seq1, kseq_t *seq2, int offset)
{
	return lex_strlt(seq1->seq.s + offset, seq2->seq.s + offset);
}


// MSEQ Utilities
/*
 * @func mask_mseq_chars
 * Masks the seq and qual fields with seqchar and qualchar.
 * Primarily used for masking barcode sequence.
 * :param: seq [mseq_t *] mseq argument.
 * :param: n_len [int] Number of characters in qual and seq to mask.
 * :param: seqchar [char] character with which to mask sequence.
 * :param: qualchar [char] character with which to mask quality.
 */
#define mask_mseq_chars(seqvar, n_len, seqchar, qualchar)\
	do {\
		memset(seqvar->seq, seqchar, n_len);\
		memset(seqvar->qual, qualchar, n_len);\
	} while(0)

/*
 * @func mask_mseq
 * Masks the seq and qual fields with 'N' and '#'.
 * :param: seq [mseq_t *] mseq argument.
 * :param: n_len [int] Number of characters in qual and seq to mask.
 */
#define mask_mseq(seqvar, n_len) mask_mseq_chars(seqvar, n_len, 'N', '#')

void mseq_destroy(mseq_t *mvar);
mseq_t *mseq_init(kseq_t *seq, char *rescaler, int is_read2);
mseq_t *mseq_rescale_init(kseq_t *seq, char *rescaler, tmp_mseq_t *tmp, int is_read2);
void update_mseq(mseq_t *mvar, kseq_t *seq, char *rescaler, tmp_mseq_t *tmp, int n_len, int is_read2, int switch_reads);
static inline void mseq2fq_stranded(FILE *handle, mseq_t *mvar, int pass_fail, char *barcode, char prefix)
{
	fprintf(handle, "@%s ~#!#~|FP=%c|BS=%c%s\n%s\n+\n%s\n",
			mvar->name, pass_fail + '0', prefix, barcode, mvar->seq, mvar->qual);
}

static inline void mseq2fq_inline(FILE *handle, mseq_t *mvar, int pass_fail, char *barcode)
{
	fprintf(handle, "@%s ~#!#~|FP=%c|BS=Z%s\n%s\n+\n%s\n",
			mvar->name, pass_fail + '0', barcode, mvar->seq, mvar->qual);
}

// TMP_MSEQ Utilities

tmp_mseq_t *init_tm_ptr(int readlen, int blen);
void tm_destroy(tmp_mseq_t *var);

#endif /* MSEQ_H */
