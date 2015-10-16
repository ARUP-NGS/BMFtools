#ifndef KSEQ_DEC_H
#define KSEQ_DEC_H

#include "kseq.h"
#include "kingfisher.h"
#include "o_mem.h"
#include <zlib.h>
/*
 * Utilities regarding kseq types.
 */
KSEQ_INIT(gzFile, gzread)

char nuc_cmpl(char character);
int nuc_cmp(char forward, char reverse);



// mseq is a mutable struct holding kseq's information.

typedef struct mseq {
	char name[100];
	char comment[2000];
	char seq[200];
	char qual[200];
	char barcode[MAX_BARCODE_LENGTH + 1];
	int l;
	int blen;
	char rc;
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
inline char *barcode_mem_view(kseq_t *seq)
{
	int hits = 0;
	for(int i = 0; i < seq->comment.l; ++i) {
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

void update_mseq(mseq_t *mvar, char *barcode, kseq_t *seq, char *rescaler, tmp_mseq_t *tmp, int n_len, int is_read2);
tmp_mseq_t *init_tm_ptr(int readlen, int blen);
void tm_destroy(tmp_mseq_t *var);
tmp_mseq_t init_tmp_mseq(int readlen, int blen);
void mseq_destroy(mseq_t *mvar);
void crc_mseq(mseq_t *mvar, tmp_mseq_t *tmp);
int crc_flip(mseq_t *mvar, char *barcode, int blen, int readlen);

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
		cmp_ret = nuc_cmp(mvar->seq[i], mvar->seq[readlen - i - 1]);
		if(cmp_ret < 0) {
			return 0; // It's lexicographically lower as is. Don't flip!
		}
		else if(cmp_ret > 0) {
			return 1; // It's not lexicographically lower as it is. Flip!
		}
	}
	return 0; // Both barcode and read are lexicographically identical forward and reverse... is that possible? Eh. Don't flip.
}

/*
 * Init's tmp_mseq_t ptr
 */
inline tmp_mseq_t *init_tm_ptr(int readlen, int blen)
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

inline void tm_destroy(tmp_mseq_t *var) {
	cond_free(var->tmp_barcode);
	cond_free(var->tmp_seq);
	cond_free(var->tmp_qual);
	cond_free(var);
}


inline tmp_mseq_t init_tmp_mseq(int readlen, int blen)
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
	cond_free(mvar.tmp_seq);
	cond_free(mvar.tmp_qual);
	cond_free(mvar.tmp_barcode);
	mvar.readlen = 0;
	mvar.blen = 0;
}

inline void mseq2fq_inline(FILE *handle, mseq_t *mvar, char pass_fail)
{
	fprintf(handle, "@%s ~#!#~|FP=%c|BS=%s|RC=%c\n%s\n+\n%s\n",
			mvar->name, pass_fail, mvar->barcode, mvar->rc, mvar->seq, mvar->qual);
	return;
}


inline void crc_mseq(mseq_t *mvar, tmp_mseq_t *tmp)
{
	if(!crc_flip(mvar, mvar->barcode, tmp->blen, tmp->readlen)) {
		mvar->rc = '0';
		return;
	}
	mvar->rc = '1';
	for(int i = 0; i < tmp->readlen; i++) {
		NUC_CMPL(mvar->seq[tmp->readlen - i - 1], tmp->tmp_seq[i])
		// Equivalent to
		//tmp->tmp_seq[i] = nuc_cmpl(mvar->seq[tmp->readlen - i - 1]);
		tmp->tmp_qual[i] = mvar->qual[tmp->readlen - i - 1];
	}
	for(int i = 0; i < tmp->blen; i++) {
		NUC_CMPL(mvar->barcode[tmp->blen - i - 1], tmp->tmp_barcode[i]);
		// Equivalent to
		//tmp->tmp_barcode[i] = nuc_cmpl(mvar->barcode[tmp->blen - i - 1]);
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


/*
 * :param: [kseq_t *] seq - kseq handle
 * :param: [mseq_t *] ret - initialized mseq_t pointer.
 * :param: [char ****] rescaler - pointer to a 3-dimensional array of rescaled phred scores.
 * :param: [tmp_mseq_t *] tmp - pointer to a tmp_mseq_t object
 * for holding information for conditional reverse complementing.
 * :param: [int] n_len - the number of bases to N at the beginning of each read.
 * :param: [int] is_read2 - true if the read is read2. Assumption: is_read2 is 0 or 1.
 */
inline mseq_t *p7_mseq_rescale_init(kseq_t *seq, char *rescaler, int n_len, int is_read2)
{
	if(!seq) {
		fprintf(stderr, "kseq for initiating p7_mseq is null. Abort!\n");
		exit(EXIT_FAILURE);
	}
	mseq_t *ret = (mseq_t *)malloc(sizeof(mseq_t));
	strcpy(ret->name, seq->name.s);
	strcpy(ret->comment, seq->comment.s);
	strcpy(ret->seq, seq->seq.s);
	//if(!rescaler) {
		strcpy(ret->qual, seq->qual.s);
	//}

	if(rescaler) {
		//ret->qual = (char *)malloc((seq->seq.l + 1) * sizeof(char));
		ret->qual[seq->seq.l] = '\0'; // Null-terminate this string.
		for(int i = 0; i < seq->seq.l; i++) {
			//fprintf(stderr, "RS params: %i, %i, %i, %i, %i.\n", is_read2 , seq->qual.s[i], i, ret->seq[i], seq->seq.l);
			ret->qual[i] = rescale_qscore(is_read2 , seq->qual.s[i], i, ret->seq[i], seq->seq.l, rescaler);
			//fprintf(stderr, "New qscore: %i.\n", ret->qual[i]);
		}
	}
	memset(ret->seq, 'N', n_len); // Set the beginning of the read to Ns.
	memset(ret->qual, '#', n_len); // Set all N bases to quality score of 2.
	ret->l = seq->seq.l;
	memset(ret->barcode, 0, MAX_BARCODE_LENGTH + 1);
	return ret;
}

/*
 * :param: [kseq_t *] seq - kseq handle
 * :param: [mseq_t *] ret - initialized mseq_t pointer.
 * :param: [char ****] rescaler - pointer to a 3-dimensional array of rescaled phred scores.
 * :param: [tmp_mseq_t *] tmp - pointer to a tmp_mseq_t object
 * for holding information for conditional reverse complementing.
 * :param: [int] n_len - the number of bases to N at the beginning of each read.
 * :param: [int] is_read2 - true if the read is read2.
 */
inline mseq_t *mseq_rescale_init(kseq_t *seq, char *rescaler, tmp_mseq_t *tmp, int n_len, int is_read2)
{
	mseq_t *ret = p7_mseq_rescale_init(seq, rescaler, n_len, is_read2);
	//fprintf(stderr, "Pointer to ret: %p. To tmp: %p. Barcode: %s.\n", ret, tmp, ret->barcode);
	if(!tmp) {
		fprintf(stderr, "Tmpvars not allocated!\n");
		exit(EXIT_FAILURE);
	}
	ret->blen = tmp->blen;
	ret->rc = 0;
	if(ret) {
		crc_mseq(ret, tmp);
	}
	return ret;
}

/*
 * Set is_read2 to 1 for read 2, 0 for read 1.
 */
inline void update_mseq(mseq_t *mvar, char *barcode, kseq_t *seq, char *rescaler, tmp_mseq_t *tmp, int n_len, int is_read2)
{
	memcpy(mvar->name, seq->name.s, seq->name.l);
    mvar->name[seq->name.l] = '\0';
	memcpy(mvar->seq, seq->seq.s, seq->seq.l * sizeof(char));
	memset(mvar->seq, 'N', n_len);
	if(!rescaler) {
		memcpy(mvar->qual, seq->qual.s, seq->qual.l * sizeof(char));
	}
	else {
		for(int i = n_len; i < seq->seq.l; i++) {
			// Leave quality scores alone for bases which are N. Otherwise
			mvar->qual[i] = (mvar->seq[i] == 'N') ? 33 : rescale_qscore(is_read2, seq->qual.s[i], i, mvar->seq[i], seq->seq.l, rescaler);
		}
	}
#if !NDEBUG
	if(strlen(mvar->qual) != seq->qual.l){
		fprintf(stderr, "Ret qual has the wrong length. (%i). Expected: %i. Seq: %s. Kseq: %s.\n", strlen(mvar->qual), seq->qual.l, mvar->qual, seq->qual.s);
		exit(1);
	}
#endif
	strcpy(mvar->barcode, barcode);
	crc_mseq(mvar, tmp);
}

inline mseq_t *init_crms_mseq(kseq_t *seq, char *barcode, char *rescaler, tmp_mseq_t *tmp, int n_len, int is_read2)
{
#if !NDEBUG
	if(!barcode) {
		fprintf(stderr, "Barocde is NULL ABORT>asdfasfjafjhaksdfkjasdfas.\n");
		exit(1);
	}
	else {
		fprintf(stderr, "Barcode: %s.\n", barcode);
	}
	if(rescaler) {
		for(int i = 0; i < seq->seq.l * 39 * 4 * 2; ++i) {
			if(rescaler[i] < 0) {
				fprintf(stderr, "Rescaler's got a negative number in init_crms_mseq. WTF? %i. Index: %i.\n", rescaler[i], i);
				exit(EXIT_FAILURE);
			}
			else if(!rescaler[i]) {
				fprintf(stderr, "Rescaler's got a zero value in init_crms_mseq. WTF? %i. Index: %i.\n", rescaler[i], i);
				exit(EXIT_FAILURE);
			}
			else {
				fprintf(stderr, "Rescaler's looking like it's supposed to. %i.Index: %i\n", rescaler[i], i);
			}
		}
	}
	fprintf(stderr, "Finished checking the array values. Now initializing mseq_t for read %i.\n", is_read2 + 1);
#endif
	mseq_t *ret = mseq_rescale_init(seq, rescaler, tmp, n_len, is_read2);
	strcpy(ret->barcode, barcode);
	return ret;
}

inline mseq_t init_rescale_revcmp_mseq(kseq_t *seq, char *barcode, char *rescaler, tmp_mseq_t *tmp, int n_len, int is_read2)
{
	mseq_t ret = {
			.l = 0,
			.blen = 0,
			.rc = '0'
	};
	memcpy(ret.barcode, barcode, strlen(barcode) + 1);
	mseq_rescale_init(seq, rescaler, tmp, n_len, is_read2);
	return ret;
}


inline void mseq_destroy(mseq_t *mvar)
{
	// Note: do not free barcode, as that is owned by another.
	mvar->l = 0;
	mvar->blen = 0;
	cond_free(mvar);
	return;
}


inline void pushback_rescaled_kseq(KingFisher_t *kfp, kseq_t *seq, char *rescaler, int *nuc_indices, int blen, int is_read2)
{
	for(int i = 0; i < kfp->readlen; i++) {
		nuc_to_pos((seq->seq.s[i]), nuc_indices);
		kfp->nuc_counts[i][nuc_indices[0]] += 1;
		kfp->phred_sums[i][nuc_indices[1]] += rescale_qscore(is_read2 ? 1 : 0, seq->qual.s[i], i, seq->seq.s[i], seq->seq.l, rescaler);
		if(seq->qual.s[i] > kfp->max_phreds[i]) {
			kfp->max_phreds[i] = seq->qual.s[i];
		}
	}
	if(kfp->length == 0) {
		char *bs_ptr = barcode_mem_view(seq);
		kfp->pass_fail = (char)*(bs_ptr- 5);
		memcpy(kfp->barcode, bs_ptr, blen);
		kfp->barcode[blen] = '\0';
	}
	kfp->length++; // Increment
	return;
}


inline void set_barcode(kseq_t *seq1, kseq_t *seq2, char *barcode, int offset, int blen1_2)
{
	memcpy(barcode, seq1->seq.s + offset, blen1_2 * sizeof(char)); // Copying the fist half of the barcode
	memcpy(barcode + blen1_2, seq2->seq.s + offset,
		   blen1_2 * sizeof(char));
	barcode[blen1_2 * 2] = '\0';
	return;
}


inline void pushback_kseq(KingFisher_t *kfp, kseq_t *seq, int *nuc_indices, int blen)
{
	for(int i = 0; i < kfp->readlen; i++) {
		nuc_to_pos((seq->seq.s[i]), nuc_indices);
		kfp->nuc_counts[i][nuc_indices[1]] += 1;
		kfp->phred_sums[i][nuc_indices[0]] += seq->qual.s[i] - 33;
		if(seq->qual.s[i] > kfp->max_phreds[i]) {
			kfp->max_phreds[i] = seq->qual.s[i];
		}
	}
	if(!kfp->length) { // Empty KingFisher
		char *bs_ptr = barcode_mem_view(seq);
		kfp->pass_fail = (char)*(bs_ptr- 5);
		memcpy(kfp->barcode, bs_ptr, blen);
		kfp->barcode[blen] = '\0';
		switch(*(bs_ptr + blen + 4)) {
			case '1': ++kfp->n_rc; break;
			case '-': kfp->n_rc = INT_MIN; break;
			default: break;
		}
	}
	++kfp->length; // Increment
	return;
}

#endif
