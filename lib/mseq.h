#ifndef MSEQ_H
#define MSEQ_H
#include <stdio.h>
#include <stdint.h>
#include <zlib.h>
#include "htslib/kseq.h"
#include "dlib/cstr_util.h"
#include "dlib/compiler_util.h"
#include "dlib/mem_util.h"
#include "lib/rescaler.h"


#define MAX_BARCODE_LENGTH 30 // Maximum expected inline barcode

// Struct definitions


struct mseq_t {
    char name[100];
    char comment[2000];
    char seq[300];
    char qual[300];
    char barcode[MAX_BARCODE_LENGTH + 1];
    int l;
    int blen;
};

typedef struct tmp_mseq {
    char *tmp_seq;
    char *tmp_qual;
    char *tmp_barcode;
    int readlen;
    int blen;
} tmp_mseq_t;

// KSEQ Utilities

#ifndef KSEQ_DEC_GZ
#define KSEQ_DEC_GZ
KSEQ_INIT(gzFile, gzread)
#endif

CONST static inline char rescale_qscore(int readnum, char qscore, int cycle, char base, int readlen, char *rescaler);
CONST static inline char *mem_view(char *comment)
{
    int hits = 0;
    for(;;) {
        switch(*comment++) {
            case '|': case '\0':
                if(hits)
                    return (char *)comment + 3;
                else hits = 1; // + 3 for |BS= minus 1, since we already incremented for the switch.
        }
    }
    return NULL; // This shouldn't ever happen.
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

/*
 * :param: kseq_t *seq1 - fastq kseq handle
 * :param: kseq_t *seq2 - fastq kseq handle
 * :param: char *barcode - buffer set by function
 * :param: int offset - number of bases to skip at the start of each read
 * :param: blen1_2 - number of bases to steal from each read
 * :returns: int - whether or not it was switched == 0.
 */
static inline int set_barcode(kseq_t *seq1, kseq_t *seq2, char *barcode, int offset, int blen1_2)
{
    if(switch_test(seq1, seq2, offset)) { // seq1's barcode is lower. No switching.
        memcpy(barcode, seq1->seq.s + offset, blen1_2 * sizeof(char)); // Copying the first half of the barcode
        memcpy(barcode + blen1_2, seq2->seq.s + offset,
                blen1_2 * sizeof(char));
        barcode[blen1_2 * 2] = '\0';
        return 0;
    } else {
        memcpy(barcode, seq2->seq.s + offset, blen1_2 * sizeof(char)); // Copying the first half of the barcode
        memcpy(barcode + blen1_2, seq1->seq.s + offset,
                blen1_2 * sizeof(char));
        barcode[blen1_2 * 2] = '\0';
        return 1;
    }
}

// calls mem_view on the comment field of the kseq_t struct.
#define barcode_mem_view(seq) mem_view(seq->comment.s)



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
static inline void mseq2fq_stranded(gzFile handle, mseq_t *mvar, int pass_fail, char *barcode, char prefix)
{
    gzprintf(handle, "@%s ~#!#~|FP=%c|BS=%c%s\n%s\n+\n%s\n",
            mvar->name, pass_fail + '0', prefix, barcode, mvar->seq, mvar->qual);
}

static inline void mseq2fq(gzFile handle, mseq_t *mvar, int pass_fail, char *barcode)
{
    gzprintf(handle, "@%s ~#!#~|FP=%c|BS=Z%s\n%s\n+\n%s\n",
            mvar->name, pass_fail + '0', barcode, mvar->seq, mvar->qual);
}


/*
 * :param: [kseq_t *] seq - kseq handle
 * :param: [mseq_t *] ret - initialized mseq_t pointer.
 * :param: [char *] rescaler - pointer to a 1-dimensional projection of a 4-dimensional array of rescaled phred scores.
 * :param: [tmp_mseq_t *] tmp - pointer to a tmp_mseq_t object
 * for holding information for conditional reverse complementing.
 * :param: [int] n_len - the number of bases to N at the beginning of each read.
 * :param: [int] is_read2 - true if the read is read2.
 */
static inline void update_mseq(mseq_t *mvar, kseq_t *seq, char *rescaler, tmp_mseq_t *tmp, int n_len, int is_read2)
{
    memcpy(mvar->name, seq->name.s, seq->name.l);
	mvar->name[seq->name.l] = '\0';
    memcpy(mvar->seq, seq->seq.s + n_len, seq->seq.l - n_len);
	mvar->seq[seq->seq.l - n_len] = '\0';
	mvar->qual[seq->qual.l - n_len] = '\0';
    if(rescaler)
        for(unsigned i = n_len; i < seq->seq.l; ++i)
            mvar->qual[i - n_len] = rescale_qscore(is_read2, seq->qual.s[i], i,
                                                   seq->seq.s[i], seq->seq.l, rescaler);
    else memcpy(mvar->qual, seq->qual.s + n_len, seq->qual.l - n_len);
}

// TMP_MSEQ Utilities

tmp_mseq_t *init_tm_ptr(int readlen, int blen);
void tm_destroy(tmp_mseq_t *var);

#endif /* MSEQ_H */
