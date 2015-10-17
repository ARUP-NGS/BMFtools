#ifndef KINGFISHER_H
#define KINGFISHER_H

#include "kseq.h"
#include "stdio.h"
#include "math.h"
#include "charcmp.h"
#include "khash.h"
#include "uthash.h"
#include "o_mem.h"
#include <zlib.h>

#ifndef MAX_BARCODE_LENGTH
#define MAX_BARCODE_LENGTH 30
#endif

KSEQ_INIT(gzFile, gzread)


// Memory costs
// max_phreds --> 1 * readlen + 8 (ptr)
// barcode --> [MAX_BARCODE_LENGTH + 1] + 8 (ptr)
// pass_fail --> 1
// nuc_counts --> readlen * 8 (ptrs) + 8 (ptr) + (readlen * 5 * sizeof(int)) [20] --> 6 * readlen + 1
// (nuc_counts = 28 * readlen + 8
// phred_sums --> readlen * 8 (ptrs) + 8 (ptr) + (readlen * 4 * sizeof(double)) [40 * readlen + 8]
// readlen --> 4
// length --> 1
// 69 * readlen + 46
// 4900 + 49 --> 5kB per barcode



typedef struct tmpbuffers {
	char name_buffer[120];
	char PVBuffer[1000];
	char FABuffer[1000];
	char cons_seq_buffer[300];
	int cons_quals[300];
	int agrees[300];
} tmpbuffers_t;

typedef struct tmpvars {
	char *bs_ptr;
	int blen;
	int readlen;
	int nuc_indices[2];
	char key[MAX_BARCODE_LENGTH + 1];
	int l; // For holding ret value for seq.
	tmpbuffers_t *buffers;
} tmpvars_t;


typedef struct KingFisher {
	uint16_t **nuc_counts; // Count of nucleotides of this form
	uint32_t **phred_sums; // Sums of -10log10(p-value)
	int length; // Number of reads in family
	int readlen; // Length of reads
	char *max_phreds; // Maximum phred score observed at position. Use this as the final sequence for the quality to maintain compatibility with GATK and other tools.
	char barcode[MAX_BARCODE_LENGTH + 1];
	char pass_fail;
	int n_rc;
} KingFisher_t;


typedef struct HashKing {
	UT_hash_handle hh;
	char id[MAX_BARCODE_LENGTH + 1];
	KingFisher_t *value;
} HashKing_t;


extern double igamc(double a, double x);
//void p7_mseq_rescale_init(kseq_t *seq, mseq_t *ret, char *rescaler, int n_len, int is_read2);
static inline char rescale_qscore(int readnum, int qscore, int cycle, char base, int readlen, char *rescaler);

//Multiply a phred score by this to convert a -10log_10(x) to a -2log_e(x)
#define LOG10E_X5_INV 0.460517018598809136803598290936872841520220297725754595206665580193514521935470496
#define LOG10E_X5_1_2 0.230258509299404568401799145468436420760110148862877297603332790096757260967735248
//such as in the following macro:
#define LOG10_TO_CHI2(x) (x) * LOG10E_X5_INV
#define AVG_LOG_TO_CHI2(x) (x) * LOG10E_X5_1_2



static inline int pvalue_to_phred(double pvalue)
{
	return (int)(-10 * log10(pvalue));
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


static inline KingFisher_t init_kf(int readlen)
{
	uint16_t **nuc_counts = (uint16_t **)malloc(readlen * sizeof(uint16_t *));
	uint32_t **phred_sums = (uint32_t **)malloc(sizeof(uint32_t *) * readlen);
	for(int i = 0; i < readlen; i++) {
		nuc_counts[i] = (uint16_t *)calloc(5, sizeof(uint16_t)); // One each for A, C, G, T, and N
		phred_sums[i] = (uint32_t *)calloc(4, sizeof(uint32_t)); // One for each nucleotide
	}
	KingFisher_t fisher = {
		.nuc_counts = nuc_counts,
		.phred_sums = phred_sums,
		.length = 0,
		.readlen = readlen,
		.max_phreds = (char *)calloc(readlen + 1, 1), // Keep track of the maximum phred score observed at position.
		.n_rc = 0,
		.pass_fail = '1'
	};
	return fisher;
}


static inline void destroy_kf(KingFisher_t *kfp)
{
	for(int i = 0; i < kfp->readlen; ++i) {
		/*
#if !NDEBUG
		fprintf(stderr, "Starting to destroy.\n");
		fprintf(stderr, "Freeing nuc_counts and phred_sums %i.", i);
#endif
		 */
		free(kfp->nuc_counts[i]);
		free(kfp->phred_sums[i]);
	}
	free(kfp->nuc_counts);
	free(kfp->phred_sums);
	free(kfp->max_phreds);
}


static inline void clear_kf(KingFisher_t *kfp)
{
	for(int i = 0; i < kfp->readlen; i++) {
		memset(kfp->nuc_counts[i], 0, 5 * sizeof(int)); // And these.
		memset(kfp->phred_sums[i], 0, 4 * sizeof(uint32_t)); // Sets these to 0.
	}
	memset(kfp->max_phreds, 0, kfp->readlen); //Turn it back into an array of nulls.
	kfp->length = 0;
	return;
}


static inline int ARRG_MAX(KingFisher_t *kfp, int index)
{
	if(kfp->phred_sums[index][3] > kfp->phred_sums[index][2] &&
	   kfp->phred_sums[index][3] > kfp->phred_sums[index][1] &&
	   kfp->phred_sums[index][3] > kfp->phred_sums[index][0]) {
		return 3;
	}
	else if(kfp->phred_sums[index][2] > kfp->phred_sums[index][1] &&
			kfp->phred_sums[index][2] > kfp->phred_sums[index][0]) {
		return 2;
	}
	else if(kfp->phred_sums[index][1] > kfp->phred_sums[index][0]) {
		return 1;
	}
	else {
		return 0;
	}
}

static inline char ARRG_MAX_TO_NUC(int argmaxret)
{
	switch (argmaxret) {
		case 1: return 'C';
		case 2: return 'G';
		case 3: return 'T';
		default: return 'A';
	}
}


static inline void fill_csv_buffer(int readlen, int *arr, char *buffer, char *prefix, char typecode)
{
	char tmpbuf[20];
	sprintf(buffer, "%s%c", prefix, typecode);
	for(int i = 0; i < readlen; i++) {
		sprintf(tmpbuf, ",%i", arr[i]);
		strcat(buffer, tmpbuf);
	}
}


static inline void fill_pv_buffer(KingFisher_t *kfp, int *phred_values, char *buffer)
{
	fill_csv_buffer(kfp->readlen, phred_values, buffer, "PV:B:", 'I');
	return;
}


static inline void fill_fa_buffer(KingFisher_t *kfp, int *agrees, char *buffer)
{
	fill_csv_buffer(kfp->readlen, agrees, buffer, "FA:B:", 'I'); // Add in the "I" to type the array.
	return;
}

/*
static inline void fill_csv_buffer_fs1(int readlen, int *arr, char *buffer, char *prefix, char typecode)
{
	char tmpbuf[20];
	sprintf(buffer, "%s%c", prefix, typecode);
	for(int i = 0; i < readlen; i++) {
		strcat(buffer, ",1");
	}
}


static inline void fill_fa_buffer_fs1(KingFisher_t *kfp, int *agrees, char *buffer)
{
	fill_csv_buffer_fs1(kfp->readlen, agrees, buffer, "FA:B:", 'I'); // Add in the "I" to type the array.
	return;
}


static inline void dmp_process_write_fs1(KingFisher_t *kfp, FILE *handle, int blen, tmpbuffers_t *tmp)
{
	//1. Argmax on the phred_sums arrays, using that to fill in the new seq and
	//buffer[0] = '@'; Set this later?
	int argmaxret;
	tmp->cons_seq_buffer[kfp->readlen] = '\0'; // Null-terminal cons_seq.
	for(int i = 0; i < kfp->readlen; ++i) {
		argmaxret = ARRG_MAX(kfp, i);
		tmp->cons_quals[i] = kfp->phred_sums[i][argmaxret];
		// Final quality must be 2 or greater and at least one read in the family should support that base call.
		tmp->cons_seq_buffer[i] = (tmp->cons_quals[i] > 2 && kfp->nuc_counts[i][argmaxret]) ? ARRG_MAX_TO_NUC(argmaxret): 'N';
		tmp->agrees[i] = kfp->nuc_counts[i][argmaxret];
	}
	fill_fa_buffer_fs2(kfp, tmp->agrees, tmp->FABuffer);
	//fprintf(stderr, "FA buffer: %s.\n", FABuffer);
	fill_pv_buffer(kfp, tmp->cons_quals, tmp->PVBuffer);
	tmp->name_buffer[0] = '@';
	memcpy((char *)(tmp->name_buffer + 1), kfp->barcode, blen);
	tmp->name_buffer[1 + blen] = '\0';
	//fprintf(stderr, "Name buffer: %s\n", tmp->name_buffer);
	//fprintf(stderr, "Output result: %s %s", tmp->name_buffer, arr_tag_buffer);
	fprintf(handle, "%s %s\t%s\tFP:i:%c\tRC:i:%i\tFM:i:%i\n%s\n+\n%s\n", tmp->name_buffer,
			tmp->FABuffer, tmp->PVBuffer,
			kfp->pass_fail, kfp->n_rc, kfp->length,
			tmp->cons_seq_buffer, kfp->max_phreds);
	return;
}
*/


/*
 * This returns primarily negative numbers. Whoops.
 */
static inline void dmp_process_write_full_pvalues(KingFisher_t *kfp, FILE *handle, int blen, tmpbuffers_t *tmp)
{
	//1. Argmax on the phred_sums arrays, using that to fill in the new seq and
	//buffer[0] = '@'; Set this later?
	int argmaxret;
	double tmp_phred;
	tmp->cons_seq_buffer[kfp->readlen] = '\0'; // Null-terminal cons_seq.
	for(int i = 0; i < kfp->readlen; ++i) {
		argmaxret = ARRG_MAX(kfp, i);
		tmp_phred = kfp->phred_sums[i][argmaxret];
		for(int j = 0; j != argmaxret && j < 4; ++j) {
			tmp_phred -= kfp->phred_sums[i][j];
		}
		tmp->cons_quals[i] = tmp_phred > 0 ? pvalue_to_phred(LOG10_TO_CHI2(tmp_phred)): 0;
		if(tmp->cons_quals[i] < -1073741824) { // Underflow!
			tmp->cons_quals[i] = 6666;
		}
		/*
		else if(tmp->cons_quals[i] < 0) {
			tmp->cons_quals[i] = 0;
		}
		*/

		// Final quality must be 2 or greater and at least one read in the family should support that base call.
		tmp->cons_seq_buffer[i] = (tmp->cons_quals[i] > 2 && kfp->nuc_counts[i][argmaxret]) ? ARRG_MAX_TO_NUC(argmaxret): 'N';
		tmp->agrees[i] = kfp->nuc_counts[i][argmaxret];
	}
	fill_fa_buffer(kfp, tmp->agrees, tmp->FABuffer);
	//fprintf(stderr, "FA buffer: %s.\n", FABuffer);
	fill_pv_buffer(kfp, tmp->cons_quals, tmp->PVBuffer);
	tmp->name_buffer[0] = '@';
	memcpy((char *)(tmp->name_buffer + 1), kfp->barcode, blen);
	tmp->name_buffer[1 + blen] = '\0';
	//fprintf(stderr, "Name buffer: %s\n", tmp->name_buffer);
	//fprintf(stderr, "Output result: %s %s", tmp->name_buffer, arr_tag_buffer);
	fprintf(handle, "%s %s\t%s\tFP:i:%c\tRC:i:%i\tFM:i:%i\n%s\n+\n%s\n", tmp->name_buffer,
			tmp->FABuffer, tmp->PVBuffer,
			kfp->pass_fail, kfp->n_rc, kfp->length,
			tmp->cons_seq_buffer, kfp->max_phreds);
	return;
}


/*
 * TODO: Use tmpvals_t object to avoid allocating and deallocating each of these.
 */
static inline void dmp_process_write_sub_chi2(KingFisher_t *kfp, FILE *handle, int blen, tmpbuffers_t *tmp)
{
	//1. Argmax on the phred_sums arrays, using that to fill in the new seq and
	//buffer[0] = '@'; Set this later?
	uint32_t tmp_max;
	int tmp_phreds[4];
	tmp->cons_seq_buffer[kfp->readlen] = '\0'; // Null-terminal cons_seq.
	for(int i = 0; i < kfp->readlen; ++i) {
		tmp_max = -1;
		for(int j = 0; j < 4; ++j) {
			if(kfp->phred_sums[i][j] > tmp_max) {
				tmp_max = j;
			}
			tmp_phreds[j] = pvalue_to_phred(igamc_pvalues(kfp->nuc_counts[i][j], LOG10_TO_CHI2(kfp->phred_sums[i][j])));
		}
		if(tmp_max < 0) {
			tmp_max = 0;
		}
		for(int j = 0; j < 4 && j != tmp_max; ++j) {
			tmp_phreds[tmp_max] -= tmp_phreds[j];
		}
		tmp->cons_quals[i] = tmp_phreds[tmp_max];
		if(tmp->cons_quals[i] < -1073741824) { // Underflow!
			tmp->cons_quals[i] = 3114;
		}
		// Final quality must be 2 or greater and at least one read in the family should support that base call.
		tmp->cons_seq_buffer[i] = (tmp->cons_quals[i] > 2 && kfp->nuc_counts[i][tmp_max]) ? ARRG_MAX_TO_NUC(tmp_max): 'N';
		tmp->agrees[i] = kfp->nuc_counts[i][tmp_max];
	}
	fill_fa_buffer(kfp, tmp->agrees, tmp->FABuffer);
	//fprintf(stderr, "FA buffer: %s.\n", FABuffer);
	fill_pv_buffer(kfp, tmp->cons_quals, tmp->PVBuffer);
	tmp->name_buffer[0] = '@';
	memcpy((char *)(tmp->name_buffer + 1), kfp->barcode, blen);
	tmp->name_buffer[1 + blen] = '\0';
	//fprintf(stderr, "Name buffer: %s\n", tmp->name_buffer);
	//fprintf(stderr, "Output result: %s %s", tmp->name_buffer, arr_tag_buffer);
	fprintf(handle, "%s %s\t%s\tFP:i:%c\tRC:i:%i\tFM:i:%i\n%s\n+\n%s\n", tmp->name_buffer,
			tmp->FABuffer, tmp->PVBuffer,
			kfp->pass_fail, kfp->n_rc, kfp->length,
			tmp->cons_seq_buffer, kfp->max_phreds);
	return;
}


/*
 * TODO: Use tmpvals_t object to avoid allocating and deallocating each of these.
 */
static inline void dmp_process_write(KingFisher_t *kfp, FILE *handle, int blen, tmpbuffers_t *tmp)
{
	//1. Argmax on the phred_sums arrays, using that to fill in the new seq and
	//buffer[0] = '@'; Set this later?
	int argmaxret;
	tmp->cons_seq_buffer[kfp->readlen] = '\0'; // Null-terminal cons_seq.
	for(int i = 0; i < kfp->readlen; ++i) {
		argmaxret = ARRG_MAX(kfp, i);
		tmp->cons_quals[i] = pvalue_to_phred(igamc_pvalues(kfp->length, LOG10_TO_CHI2((kfp->phred_sums[i][argmaxret]))));
		if(tmp->cons_quals[i] < -1073741824) { // Underflow!
			tmp->cons_quals[i] = 3114;
		}
		// Final quality must be 2 or greater and at least one read in the family should support that base call.
		tmp->cons_seq_buffer[i] = (tmp->cons_quals[i] > 2 && kfp->nuc_counts[i][argmaxret]) ? ARRG_MAX_TO_NUC(argmaxret): 'N';
		tmp->agrees[i] = kfp->nuc_counts[i][argmaxret];
	}
	fill_fa_buffer(kfp, tmp->agrees, tmp->FABuffer);
	//fprintf(stderr, "FA buffer: %s.\n", FABuffer);
	fill_pv_buffer(kfp, tmp->cons_quals, tmp->PVBuffer);
	tmp->name_buffer[0] = '@';
	memcpy((char *)(tmp->name_buffer + 1), kfp->barcode, blen);
	tmp->name_buffer[1 + blen] = '\0';
	//fprintf(stderr, "Name buffer: %s\n", tmp->name_buffer);
	//fprintf(stderr, "Output result: %s %s", tmp->name_buffer, arr_tag_buffer);
	fprintf(handle, "%s %s\t%s\tFP:i:%c\tRC:i:%i\tFM:i:%i\n%s\n+\n%s\n", tmp->name_buffer,
			tmp->FABuffer, tmp->PVBuffer,
			kfp->pass_fail, kfp->n_rc, kfp->length,
			tmp->cons_seq_buffer, kfp->max_phreds);
	return;
}


static inline char rescale_qscore(int readnum, int qscore, int cycle, char base, int readlen, char *rescaler)
{
	int index = readnum;
	int mult = 2;
	//fprintf(stderr, "index value is now: %i, mult %i.\n", index, mult);
	index += cycle * mult;
	mult *= readlen;
	//fprintf(stderr, "index value is now: %i, mult %i. Qscore: %i, Qscore index%i.\n", index, mult, qscore, qscore - 35);
	index += (qscore - 35) * mult; // Subtract 35 - 33 to get to phred space, 2 to offset by 2.
	mult *= 39;
	//fprintf(stderr, "index value is now: %i, mult %i.\n", index, mult);
	index += mult * nuc2num(base);
	//fprintf(stderr, "Index = %i.\n", index);
#if !NDEBUG
	if(index >= readlen * 2 * 39 * 4 || index < 0) {
		fprintf(stderr, "Something's wrong. Index (%i) is too big or negative! Max: %i.\n", index, readlen * 2 * 39 * 4);
		exit(EXIT_FAILURE);
	}
	//fprintf(stderr, "Value at index: %i (%c).\n", rescaler[index], rescaler[index] + 33);
	if(rescaler[index] < 0) {
		fprintf(stderr, "WTF THIS CAN'T BE BELOW 0 (%i).\n", rescaler[index]);
		exit(EXIT_FAILURE);
	}
#endif
	return rescaler[index] + 33;
}



static inline void set_kf(int readlen, KingFisher_t ret)
{
	ret.length = 0;
	ret.readlen = readlen;
	ret.nuc_counts = (uint16_t **)malloc(readlen * sizeof(uint16_t *));
	ret.phred_sums = (uint32_t **)malloc(sizeof(uint32_t *) * readlen);
	ret.max_phreds = (char *)calloc(readlen + 1, sizeof(char)); // Keep track of the maximum phred score observed at position.
	for(int i = 0; i < readlen; i++) {
		ret.nuc_counts[i] = (uint16_t *)calloc(5, sizeof(uint16_t)); // One each for A, C, G, T, and N
		ret.phred_sums[i] = (uint32_t *)calloc(4, sizeof(uint32_t)); // One for each nucleotide
	}
	return;
}


static inline KingFisher_t *init_kfp(size_t readlen)
{
	KingFisher_t *ret = (KingFisher_t *)malloc(sizeof(KingFisher_t));
	ret->length = 0; // Check to see if this is necessary after calloc - I'm pretty sure not.
	ret->n_rc = 0;
	ret->readlen = readlen;
	ret->max_phreds = (char *)calloc(readlen + 1, sizeof(char)); // Keep track of the maximum phred score observed at position.
	ret->nuc_counts = (uint16_t **)malloc(readlen * sizeof(uint16_t *));
	ret->phred_sums = (uint32_t **)malloc(readlen * sizeof(uint32_t *));
	for(int i = 0; i < readlen; ++i) {
		ret->nuc_counts[i] = (uint16_t *)calloc(5, sizeof(uint16_t)); // One each for A, C, G, T, and N
		ret->phred_sums[i] = (uint32_t *)calloc(4, sizeof(uint32_t)); // One for each nucleotide
	}
	ret->pass_fail = '1';
	return ret;
}


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
static inline char *barcode_mem_view(kseq_t *seq)
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


static inline int crc_flip(mseq_t *mvar, char *barcode, int blen, int readlen)
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
static inline tmp_mseq_t *init_tm_ptr(int readlen, int blen)
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

static inline void tm_destroy(tmp_mseq_t *var) {
	cond_free(var->tmp_barcode);
	cond_free(var->tmp_seq);
	cond_free(var->tmp_qual);
	cond_free(var);
}


static inline tmp_mseq_t init_tmp_mseq(int readlen, int blen)
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


static inline void tmp_mseq_destroy(tmp_mseq_t mvar)
{
	cond_free(mvar.tmp_seq);
	cond_free(mvar.tmp_qual);
	cond_free(mvar.tmp_barcode);
	mvar.readlen = 0;
	mvar.blen = 0;
}

static inline void mseq2fq_inline(FILE *handle, mseq_t *mvar, char pass_fail)
{
	fprintf(handle, "@%s ~#!#~|FP=%c|BS=%s|RC=%c\n%s\n+\n%s\n",
			mvar->name, pass_fail, mvar->barcode, mvar->rc, mvar->seq, mvar->qual);
	return;
}


static inline void crc_mseq(mseq_t *mvar, tmp_mseq_t *tmp)
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
static inline mseq_t *p7_mseq_rescale_init(kseq_t *seq, char *rescaler, int n_len, int is_read2)
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
static inline mseq_t *mseq_rescale_init(kseq_t *seq, char *rescaler, tmp_mseq_t *tmp, int n_len, int is_read2)
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
static inline void update_mseq(mseq_t *mvar, char *barcode, kseq_t *seq, char *rescaler, tmp_mseq_t *tmp, int n_len, int is_read2)
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

static inline mseq_t *init_crms_mseq(kseq_t *seq, char *barcode, char *rescaler, tmp_mseq_t *tmp, int n_len, int is_read2)
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

static inline mseq_t init_rescale_revcmp_mseq(kseq_t *seq, char *barcode, char *rescaler, tmp_mseq_t *tmp, int n_len, int is_read2)
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


static inline void mseq_destroy(mseq_t *mvar)
{
	// Note: do not free barcode, as that is owned by another.
	mvar->l = 0;
	mvar->blen = 0;
	cond_free(mvar);
	return;
}


static inline void pushback_rescaled_kseq(KingFisher_t *kfp, kseq_t *seq, char *rescaler, int *nuc_indices, int blen, int is_read2)
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


static inline void set_barcode(kseq_t *seq1, kseq_t *seq2, char *barcode, int offset, int blen1_2)
{
	memcpy(barcode, seq1->seq.s + offset, blen1_2 * sizeof(char)); // Copying the fist half of the barcode
	memcpy(barcode + blen1_2, seq2->seq.s + offset,
		   blen1_2 * sizeof(char));
	barcode[blen1_2 * 2] = '\0';
	return;
}


static inline void pushback_kseq(KingFisher_t *kfp, kseq_t *seq, int *nuc_indices, int blen)
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
