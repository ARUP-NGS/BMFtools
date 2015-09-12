#pragma once
#include "kseq.h"
#include "kingfisher.h"
#include <zlib.h>
/*
 * Utilities regarding kseq types.
 */
KSEQ_INIT(gzFile, gzread)

char nuc_cmpl(char character);
int nuc_cmp(char forward, char reverse);
int rescale_qscore(int qscore, int cycle, char base, char ***rescaler);



// mseq is a mutable struct holding kseq's information.

typedef struct mseq {
    char *name;
    char *comment;
    char *seq;
    char *qual;
    char *barcode;
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


/*
 * Warning: returns a NULL upon not finding a second pipe symbol.
 * This is *NOT* a properly null-terminated string.
 */
char *barcode_mem_view(kseq_t *seq)
{
    int hits = 0;
    for(int i = 0; i < seq->comment.l; i++) {
        if(seq->comment.s[i] == '|') {
            if(!hits) {
                hits += 1;
            }
            else {
                return (char *)(seq->comment.s + i + 4); // 4 for "|BS="
            }
        }
    }
    return NULL;
}


inline int crc_flip(mseq_t *mvar, char *barcode, int blen, int readlen)
{
    int cmp_ret;
    for(int i = 0; i < blen; ++i) {
        cmp_ret = nuc_cmp(barcode[i], barcode[blen - i - 1]);
        if(cmp_ret < 0) {
            return 0; // It's lexicographically lower as is. Don't flip!
        }
        else if(cmp_ret > 0) {
            return 1; // It's not lexicographically lower as it is. Flip!
        }
    }
    for(int i = 0; i < readlen; ++i) {
        cmp_ret = nuc_cmp(barcode[i], barcode[readlen - i - 1]);
        if(cmp_ret < 0) {
            return 0; // It's lexicographically lower as is. Don't flip!
        }
        else if(cmp_ret > 0) {
            return 1; // It's not lexicographically lower as it is. Flip!
        }
    }
    return 0; // Both barcode and read are lexicographically identical forward and reverse... is that possible? Eh. Don't flip.
}


tmp_mseq_t init_tmp_mseq(int readlen, int blen)
{
    char *tmp_seq = (char *)malloc(readlen * sizeof(char));
    char *tmp_qual = (char *)malloc(readlen * sizeof(char));
    char *tmp_barcode = (char *)malloc(blen * sizeof(char));
    tmp_mseq_t ret = {
        .tmp_seq = tmp_seq,
        .tmp_qual = tmp_qual,
        .tmp_barcode = tmp_barcode,
        .readlen = readlen,
        .blen = blen
    };
    return ret;
}


void tmp_mseq_destroy(tmp_mseq_t mvar)
{
    free(mvar.tmp_seq);
    free(mvar.tmp_qual);
    free(mvar.tmp_barcode);
    mvar.readlen = 0;
    mvar.blen = 0;
}

inline void mseq2fq_inline(FILE *handle, mseq_t *mvar, char pass_fail)
{
    fprintf(handle, "@%s ~#!#~|FP=%c|BS=%s\n%s\n+\n%s\n",
            mvar->name, pass_fail, mvar->barcode, mvar->seq, mvar->qual);
    return;
}


inline void crc_mseq(mseq_t *mvar, tmp_mseq_t *tmp)
{
    if(!crc_flip(mvar, mvar->barcode, tmp->blen, tmp->readlen)) return;
    for(int i = 0; i < tmp->readlen; i++) {
        tmp->tmp_seq[i] = nuc_cmpl(mvar->seq[tmp->readlen - i - 1]);
        tmp->tmp_qual[i] = mvar->qual[tmp->readlen - i - 1];
    }
    for(int i = 0; i < tmp->blen; i++) {
        tmp->tmp_barcode[i] = nuc_cmpl(mvar->barcode[tmp->blen - i - 1]);
    }
#if !NDEBUG
    char *omgzwtf = (char *)malloc(tmp->readlen + 1);
    omgzwtf[tmp->readlen] = '\0';
    memcpy(omgzwtf, tmp->tmp_seq, tmp->readlen);
    free(omgzwtf);
    omgzwtf = (char *)malloc(tmp->blen + 1);
    omgzwtf[tmp->blen] = '\0';
    memcpy(omgzwtf, tmp->tmp_barcode, tmp->blen);
    free(omgzwtf);
#endif
    memcpy(mvar->qual, tmp->tmp_qual, tmp->readlen * sizeof(char));
    memcpy(mvar->seq, tmp->tmp_seq, tmp->readlen * sizeof(char));
    memcpy(mvar->barcode, tmp->tmp_barcode, tmp->blen * sizeof(char));
    return;
}


inline void mseq_rescale_init(kseq_t *seq, mseq_t *ret, char ***rescaler, tmp_mseq_t *tmp, int n_len)
{
    if(!seq) {
        ret = NULL;
        return;
    }
    ret->name = strdup(seq->name.s);
    ret->comment = strdup(seq->comment.s);
    ret->seq = strdup(seq->seq.s);
    memset(ret->seq, 78, n_len); // Set the beginning of the read to Ns.
    if(!rescaler) ret->qual = strdup(seq->qual.s);
    else {
        ret->qual = (char *)malloc((seq->seq.l + 1) * sizeof(char));
        ret->qual[seq->seq.l] = '\0'; // Null-terminate this string.
        for(int i = 0; i < seq->seq.l; i++) {
            ret->qual[i] = rescale_qscore(ret->qual[i], i, ret->seq[i], rescaler);
        }
    }
    ret->l = seq->seq.l;
    crc_mseq(ret, tmp);
}


inline void update_mseq(mseq_t *mvar, char *barcode, kseq_t *seq, char ***rescaler, tmp_mseq_t *tmp, int n_len)
{
    memcpy(mvar->name, seq->name.s, seq->name.l * sizeof(char)); // Update name
    memcpy(mvar->seq, seq->seq.s, seq->seq.l * sizeof(char));
    memset(mvar->seq, 78, n_len);
    if(!rescaler) {
        memcpy(mvar->qual, seq->qual.s, seq->qual.l * sizeof(char));
    }
    else {
        for(int i = n_len; i < seq->seq.l; i++) {
            // Leave quality scores alone for bases which are N. Otherwise
            mvar->qual[i] = (mvar->seq[i] == 'N') ? 2 : rescale_qscore(mvar->qual[i], i, mvar->seq[i], rescaler);
        }
    }
    mvar->barcode = barcode;
    crc_mseq(mvar, tmp);
}

inline mseq_t init_rescale_revcmp_mseq(kseq_t *seq, char *barcode, char ***rescaler, tmp_mseq_t *tmp, int n_len)
{
    mseq_t ret = {
            .name = NULL,
            .comment = NULL,
            .seq = NULL,
            .qual = NULL,
            .barcode = barcode, // barcode still belongs to the argument variable!
            .l = 0,
            .blen = 0
    };
    mseq_rescale_init(seq, &ret, rescaler, tmp, n_len);
    return ret;
}


inline void mseq_destroy(mseq_t *mvar)
{
    free(mvar->name);
    free(mvar->comment);
    free(mvar->seq);
    free(mvar->qual);
    // Note: do not free barcode, as that is owned by another.
    mvar->l = 0;
    mvar->blen = 0;
    return;
}


inline void pushback_rescaled_kseq(KingFisher_t *kfp, kseq_t *seq, char ***rescaler, int *nuc_indices, int blen)
{
    fprintf(stderr, "Pushing back kseq with read length %i\n", kfp->readlen);
    for(int i = 0; i < kfp->readlen; i++) {
        nuc_to_pos((seq->seq.s[i]), nuc_indices);
        kfp->nuc_counts[i][nuc_indices[0]] += 1;
        kfp->phred_sums[i][nuc_indices[1]] += rescale_qscore(seq->qual.s[i] - 33, i, seq->seq.s[i], rescaler);
        if(seq->qual.s[i] > kfp->max_phreds[i]) {
            kfp->max_phreds[i] = seq->qual.s[i];
        }
    }
    if(kfp->length == 0) {
        char *bs_ptr = barcode_mem_view(seq);
        kfp->pass_fail = (char)*(bs_ptr- 5);
        kfp->barcode = (char *)calloc(blen + 1, sizeof(char));
        memcpy(kfp->barcode, bs_ptr, blen);
    }
    kfp->length++; // Increment
    fprintf(stderr, "New length of kfp: %i. BTW, readlen for kfp: %i.\n", kfp->length, kfp->readlen);
    return;
}


inline void set_barcode(kseq_t *seq1, kseq_t *seq2, char *barcode, int offset, int blen1_2)
{
    memcpy(barcode, seq1->seq.s + offset, blen1_2 * sizeof(char)); // Copying the fist half of the barcode
    memcpy(barcode + blen1_2, seq2->seq.s + offset,
           blen1_2 * sizeof(char));
    char *omgzwtf = (char *)malloc((blen1_2 * 2 + 1) * sizeof(char));
    omgzwtf[blen1_2 * 2] = '\0';
    memcpy(omgzwtf, barcode, blen1_2 * 2);
    free(omgzwtf);
    return;
}


inline void pushback_kseq(KingFisher_t *kfp, kseq_t *seq, int *nuc_indices, int blen)
{
    fprintf(stderr, "Pushing back kseq with read length %i\n", kfp->readlen);
    for(int i = 0; i < kfp->readlen; i++) {
        nuc_to_pos((seq->seq.s[i]), nuc_indices);
        kfp->nuc_counts[i][nuc_indices[0]] += 1;
        kfp->phred_sums[i][nuc_indices[1]] += seq->qual.s[i] - 33;
        if(seq->qual.s[i] > kfp->max_phreds[i]) {
            kfp->max_phreds[i] = seq->qual.s[i];
        }
    }
    if(kfp->length == 0) {
        char *bs_ptr = barcode_mem_view(seq);
        kfp->pass_fail = (char)*(bs_ptr- 5);
        kfp->barcode = (char *)calloc(blen + 1, sizeof(char));
        memcpy(kfp->barcode, bs_ptr, blen);
    }
    kfp->length++; // Increment
    fprintf(stderr, "New length of kfp: %i. BTW, readlen for kfp: %i.\n", kfp->length, kfp->readlen);
    return;
}

