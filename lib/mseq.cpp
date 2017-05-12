#include "mseq.h"

#include "dlib/misc_util.h"

namespace bmf {

void mseq_destroy(mseq_t *mvar)
{
    // Note: does not free barcode, as that is owned by another.
    mvar->l = 0;
    mvar->blen = 0;
    free(mvar);
}


/*
 * :param: [kseq_t *] seq - kseq handle
 * :param: [mseq_t *] ret - initialized mseq_t pointer.
 * :param: [char ****] rescaler - pointer to a 3-dimensional array of rescaled phred scores.
 * :param: [tmp_mseq_t *] tmp - pointer to a tmp_mseq_t object
 * for holding information for conditional reverse complementing.
 * :param: [int] n_len - the number of bases to N at the beginning of each read.
 * :param: [int] is_read2 - true if the read is read2. Assumption: is_read2 is 0 or 1.
 */
mseq_t *mseq_init(kseq_t *seq, char *rescaler, int is_read2)
{
    if(!seq) {
        fprintf(stderr, "kseq for initiating p7_mseq is null. Abort!\n");
        exit(EXIT_FAILURE);
    }
    mseq_t *ret((mseq_t *)calloc(1, sizeof(mseq_t)));
    strcpy(ret->name, seq->name.s);
    strcpy(ret->comment, seq->comment.s ? seq->comment.s: "");
    strcpy(ret->seq, seq->seq.s);
    strcpy(ret->qual, seq->qual.s);

    ret->l = seq->seq.l;
    if(rescaler)
        for(int i = 0; i < ret->l; i++)
            ret->qual[i] = rescale_qscore(is_read2 , seq->qual.s[i], i, ret->seq[i], seq->seq.l, rescaler);
    return ret;
}

/*
 * :param: [kseq_t *] seq - kseq handle
 * :param: [char *] rescaler - pointer to a 1-dimensional projection of a 4-dimensional array of rescaled phred scores.
 * :param: [tmp_mseq_t *] tmp - pointer to a tmp_mseq_t object
 * for holding information for conditional reverse complementing.
 * :param: [int] is_read2 - true if the read is read2.
 */
mseq_t *mseq_rescale_init(kseq_t *seq, char *rescaler, tmp_mseq_t *tmp, int is_read2)
{
    mseq_t *ret(mseq_init(seq, rescaler, is_read2));
    //fprintf(stderr, "Pointer to ret: %p. To tmp: %p. Barcode: %s.\n", ret, tmp, ret->barcode);
    if(!tmp) {
        fprintf(stderr, "Tmpvars not allocated!\n");
        exit(EXIT_FAILURE);
    }
    ret->blen = tmp->blen;
    return ret;
}


tmp_mseq_t *init_tm_ptr(int readlen, int blen)
{
    tmp_mseq_t *ret((tmp_mseq_t *)malloc(sizeof(tmp_mseq_t)));
    ret->tmp_seq = (char *)malloc(readlen * sizeof(char));
    ret->tmp_qual = (char *)malloc(readlen * sizeof(char));
    ret->tmp_barcode = (char *)malloc(blen * sizeof(char));
    ret->readlen = readlen;
    ret->blen = blen;
    return ret;
}

void tm_destroy(tmp_mseq_t *var) {
    cond_free(var->tmp_barcode);
    cond_free(var->tmp_seq);
    cond_free(var->tmp_qual);
    cond_free(var);
}

} /* namespace bmf */
