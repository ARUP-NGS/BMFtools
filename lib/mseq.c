#include "mseq.h"

CONST inline char *mem_view(char *comment)
{
	int hits = 0;
	for(;;++comment) {
		switch(*comment) == '|' || *comment == '\0') {
			case '|':
			case '\0': if(hits) return (char *)comment + 4; else hits = 1;
		}
	}
	return NULL; // This shouldn't ever happen.
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
	mseq_t *ret = (mseq_t *)calloc(1, sizeof(mseq_t));
	strcpy(ret->name, seq->name.s);
	strcpy(ret->comment, seq->comment.s);
	strcpy(ret->seq, seq->seq.s);
	strcpy(ret->qual, seq->qual.s);

	ret->l = seq->seq.l;
	if(rescaler) {
		for(int i = 0; i < ret->l; i++)
			ret->qual[i] = rescale_qscore(is_read2 , seq->qual.s[i], i, ret->seq[i], seq->seq.l, rescaler);
	}
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
	mseq_t *ret = mseq_init(seq, rescaler, is_read2);
	//fprintf(stderr, "Pointer to ret: %p. To tmp: %p. Barcode: %s.\n", ret, tmp, ret->barcode);
	if(!tmp) {
		fprintf(stderr, "Tmpvars not allocated!\n");
		exit(EXIT_FAILURE);
	}
	ret->blen = tmp->blen;
	return ret;
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
inline void update_mseq(mseq_t *mvar, kseq_t *seq, char *rescaler, tmp_mseq_t *tmp, int n_len, int is_read2, int switch_reads)
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



tmp_mseq_t *init_tm_ptr(int readlen, int blen)
{
	tmp_mseq_t *ret = (tmp_mseq_t *)malloc(sizeof(tmp_mseq_t));
	char *tmp_seq = (char *)malloc(readlen * sizeof(char));
	char *tmp_qual = (char *)malloc(readlen * sizeof(char));
	char *tmp_barcode = (char *)malloc(blen * sizeof(char));
	ret->tmp_seq = tmp_seq;
	ret->tmp_qual = tmp_qual;
	ret->tmp_barcode = tmp_barcode;
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

