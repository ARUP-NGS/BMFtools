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
/*
 * @func mem_view
 * Goes until it finds a second delimiter character, then returns a pointer 4 beyond it.
 * This points to the barcode sequence in a fastq comment.
 * :param: comment [char *] comment string
 * :returns: [char *] Pointer to the start of the barcode.
 * This is *NOT* a properly null-terminated string.
 */
CONST inline char *mem_view(char *comment)
{
	int hits = 0;
	for(;;++comment) {
		if(*comment == '|' || *comment == '\0') {
			if(hits) return (char *)comment + 4; // 4 for "|BS="
			else hits = 1;
		}
	}
	return NULL; // This shouldn't ever happen.
}
/*
 * @func mem_view
 * Goes until it finds a second delimiter character, then returns a pointer 4 beyond it.
 * This points to the barcode sequence in a fastq comment.
 * :param: seq [kseq_t *] kseq record.
 * :returns: [char *] Pointer to the start of the barcode.
 * This is *NOT* a properly null-terminated string.
 */
CONST inline char *barcode_mem_view(kseq_t *seq)
{
	return mem_view(seq->comment.s);
}

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

void mseq_destroy(mseq_t *mvar);
mseq_t *mseq_init(kseq_t *seq, char *rescaler, int is_read2);
mseq_t *mseq_rescale_init(kseq_t *seq, char *rescaler, tmp_mseq_t *tmp, int is_read2);
static inline void mseq2fq_stranded(FILE *handle, mseq_t *mvar, char pass_fail, char *barcode, char prefix)
{
	fprintf(handle, "@%s ~#!#~|FP=%c|BS=%c%s\n%s\n+\n%s\n",
			mvar->name, pass_fail, prefix, barcode, mvar->seq, mvar->qual);
}

static inline void mseq2fq_inline(FILE *handle, mseq_t *mvar, char pass_fail, char *barcode)
{
	fprintf(handle, "@%s ~#!#~|FP=%c|BS=Z%s\n%s\n+\n%s\n",
			mvar->name, pass_fail, barcode, mvar->seq, mvar->qual);
}
/*
 * :param: [kseq_t *] seq - kseq handle
 * :param: [mseq_t *] ret - initialized mseq_t pointer.
 * :param: [char *] rescaler - pointer to a 1-dimensional projection of a 4-dimensional array of rescaled phred scores.
 * :param: [tmp_mseq_t *] tmp - pointer to a tmp_mseq_t object
 * for holding information for conditional reverse complementing.
 * :param: [int] n_len - the number of bases to N at the beginning of each read.
 * :param: [int] is_read2 - true if the read is read2.
 * :param: [int] switch_reads - Whether or not to switch reads 1 and 2.
 */
static inline void update_mseq(mseq_t *mvar, kseq_t *seq, char *rescaler, tmp_mseq_t *tmp, int n_len, int is_read2, int switch_reads)
{
	memcpy(mvar->name, seq->name.s, seq->name.l);
    mvar->name[seq->name.l] = '\0';
	memcpy(mvar->seq, seq->seq.s, seq->seq.l * sizeof(char));
	memset(mvar->seq, 'N', n_len), memset(mvar->qual, '#', n_len);
	if(rescaler)
		for(int i = n_len; i < seq->seq.l; ++i)
			mvar->qual[i] = (mvar->seq[i] == 'N') ? '#' : rescale_qscore(is_read2, seq->qual.s[i], i, mvar->seq[i], seq->seq.l, rescaler);
	else
		memcpy(mvar->qual + n_len, seq->qual.s + n_len, seq->qual.l * sizeof(char) - n_len);
}


// TMP_MSEQ Utilities

tmp_mseq_t *init_tm_ptr(int readlen, int blen);
void tm_destroy(tmp_mseq_t *var);

#endif /* MSEQ_H */
