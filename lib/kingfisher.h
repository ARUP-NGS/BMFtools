#ifndef KINGFISHER_H
#define KINGFISHER_H

#include "kseq.h"
#include "stdio.h"
#include "math.h"
#include "charcmp.h"
#include "cstr_utils.h"
#include "khash.h"
#include "mem_util.h"
#include <zlib.h>
#include <inttypes.h>

#ifndef MAX_BARCODE_LENGTH
#define MAX_BARCODE_LENGTH 30
#endif

#ifndef KSEQ_DEC_GZ
#define KSEQ_DEC_GZ
KSEQ_INIT(gzFile, gzread)
#endif

#define HOM_SEQ_OFFSET 5


typedef struct tmpbuffers {
	char name_buffer[120];
	char PVBuffer[1000];
	char FABuffer[1000];
	char cons_seq_buffer[300];
	uint32_t cons_quals[300];
	uint16_t agrees[300];
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
	uint16_t *nuc_counts; // Count of nucleotides of this form
	uint32_t *phred_sums; // Sums of -10log10(p-value)
	int length; // Number of reads in family
	int readlen; // Length of reads
	char *max_phreds; // Maximum phred score observed at position. Use this as the final sequence for the quality to maintain compatibility with GATK and other tools.
	char barcode[MAX_BARCODE_LENGTH + 1];
	char pass_fail;
	int n_rc;
} KingFisher_t;


extern double igamc(double a, double x);
//void p7_mseq_rescale_init(kseq_t *seq, mseq_t *ret, char *rescaler, int n_len, int is_read2);
static inline char rescale_qscore(int readnum, int qscore, int cycle, char base, int readlen, char *rescaler);

//Multiply a phred score by this to convert a -10log_10(x) to a -2log_e(x)
#define LOG10E_X5_INV 0.460517018598809136803598290936872841520220297725754595206665580193514521935470496
#define LOG10E_X5_1_2 0.230258509299404568401799145468436420760110148862877297603332790096757260967735248
//such as in the following macro:
#define LOG10_TO_CHI2(x) (x) * LOG10E_X5_INV
#define AVG_LOG_TO_CHI2(x) (x) * LOG10E_X5_1_2



static inline uint32_t pvalue_to_phred(double pvalue)
{
	return (uint32_t)(-10 * log10(pvalue) + 0.5); // 0.5 to round up
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
	KingFisher_t fisher = {
		.nuc_counts = (uint16_t *)calloc(readlen * 5, sizeof(uint16_t)),
		.phred_sums = (uint32_t *)calloc(readlen * 4, sizeof(uint32_t)),
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
	cond_free(kfp->nuc_counts);
	cond_free(kfp->phred_sums);
	cond_free(kfp->max_phreds);
	free(kfp);
	kfp = NULL;
}


static inline void clear_kf(KingFisher_t *kfp)
{
	memset(kfp->nuc_counts, 0, 5 * sizeof(uint16_t) * kfp->readlen);
	memset(kfp->phred_sums, 0, 4 * sizeof(uint32_t) * kfp->readlen);
	memset(kfp->max_phreds, 0, kfp->readlen); //Turn it back into an array of nulls.
	kfp->length = 0;
	return;
}


static inline int ARRG_MAX(KingFisher_t *kfp, int index)
{
/*
	uint32_t max_index = 0, i;
	for(i = 0;i < 4; ++i) {
		if(kfp->phred_sums[index * 4 + i] > kfp->phred_sums[index * 4 + max_index]) {
			max_index = i;
		}
	}
	return max_index;
*/
	int i4 = index * 4;
	if(kfp->phred_sums[i4] > kfp->phred_sums[i4 + 1] &&
		kfp->phred_sums[i4] > kfp->phred_sums[i4 + 2] &&
		kfp->phred_sums[i4] > kfp->phred_sums[i4 + 3])
		return 0;

	else if(kfp->phred_sums[i4 + 1] > kfp->phred_sums[i4 + 2] &&
			kfp->phred_sums[i4 + 1] > kfp->phred_sums[i4 + 3])
		return 1;
	else if(kfp->phred_sums[i4 + 2] > kfp->phred_sums[i4 + 3])
		return 2;
	else
		return 3;
	return 0; // This never happens
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


static inline void fill_csv_buffer(int readlen, uint32_t *arr, char *buffer, char *prefix, char typecode)
{
	char tmpbuf[20];
	sprintf(buffer, "%s%c", prefix, typecode);
	for(uint32_t i = 0; i < readlen; i++) {
		sprintf(tmpbuf, ",%i", arr[i]);
		strcat(buffer, tmpbuf);
	}
}


static inline void fill_pv_buffer(KingFisher_t *kfp, uint32_t *phred_values, char *buffer)
{
	fill_csv_buffer(kfp->readlen, phred_values, buffer, "PV:B:", 'I');
	return;
}

static inline void fill_fa_buffer(KingFisher_t *kfp, uint16_t *agrees, char *buffer)
{
	char tmpbuf[20];
	sprintf(buffer, "FA:B:I");
	for(uint32_t i = 0; i < kfp->readlen; ++i) {
		sprintf(tmpbuf, ",%"PRIu16"", agrees[i]);
		strcat(buffer, tmpbuf);
	}
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



static inline void kh_pw(HashKing_t *hkp, FILE *handle, int blen, tmpbuffers_t *tmp, char *barcode)
{
	tmp->cons_seq_buffer[hkp->value->readlen] = '\0'; // Null-terminal cons_seq.
	for(int i = 0; i < hkp->value->readlen; ++i) {
		int argmaxret = ARRG_MAX(hkp->value, i);
		tmp->cons_quals[i] = pvalue_to_phred(igamc_pvalues(hkp->value->length, LOG10_TO_CHI2((hkp->value->phred_sums[i * 4 + argmaxret]))));
		// Final quality must be 2 or greater and at least one read in the family should support that base call.
		tmp->cons_seq_buffer[i] = (tmp->cons_quals[i] > 2 && hkp->value->nuc_counts[i * 4 + argmaxret]) ? ARRG_MAX_TO_NUC(argmaxret): 'N';
		tmp->agrees[i] = hkp->value->nuc_counts[i * 4 + argmaxret];
	}
	fill_fa_buffer(hkp->value, tmp->agrees, tmp->FABuffer);
	//fprintf(stderr, "FA buffer: %s.\n", FABuffer);
	fill_pv_buffer(hkp->value, tmp->cons_quals, tmp->PVBuffer);
	fprintf(handle, "@%s %s\t%s\tFP:i:%c\tRV:i:%i\tFM:i:%i\n%s\n+\n%s\n", barcode,
			tmp->FABuffer, tmp->PVBuffer,
			hkp->value->pass_fail, hkp->value->n_rc, hkp->value->length,
			tmp->cons_seq_buffer, hkp->value->max_phreds);
}
*/


static inline void dmp_process_write(KingFisher_t *kfp, FILE *handle, tmpvars_t *tmp)
{;
	//1. Argmax on the phred_sums arrays, using that to fill in the new seq and
	//buffer[0] = '@'; Set this later?
	int argmaxret;
	tmpbuffers_t *bufs = tmp->buffers;
	bufs->cons_seq_buffer[kfp->readlen] = '\0'; // Null-terminal cons_seq.
	for(int i = 0; i < kfp->readlen; ++i) {
		argmaxret = ARRG_MAX(kfp, i);
		bufs->cons_quals[i] = pvalue_to_phred(igamc_pvalues(kfp->length, LOG10_TO_CHI2((kfp->phred_sums[i * 4 + argmaxret]))));
		/*
		if(bufs->cons_quals[i] > 1073741824) { // Underflow!
			fprintf(stderr, "Note: phred_sums in underflow: %" PRIu32 ".\n", kfp->phred_sums[i * 4 + argmaxret]);
			bufs->cons_quals[i] = 3114;
		}
		*/
		// Final quality must be 2 or greater and at least one read in the family should support that base call.
		if(bufs->cons_quals[i] > 2) {
			bufs->cons_seq_buffer[i] = ARRG_MAX_TO_NUC(argmaxret);
			bufs->agrees[i] = kfp->nuc_counts[i * 4 + argmaxret];
		}
		bufs->cons_seq_buffer[i] = (bufs->cons_quals[i] > 2 && kfp->nuc_counts[i * 4 + argmaxret]) ? ARRG_MAX_TO_NUC(argmaxret): 'N';
		bufs->agrees[i] = kfp->nuc_counts[i * 4 + argmaxret];
	}
	fill_fa_buffer(kfp, bufs->agrees, bufs->FABuffer);
	//fprintf(stderr, "FA buffer: %s.\n", FABuffer);
	fill_pv_buffer(kfp, bufs->cons_quals, bufs->PVBuffer);
	fprintf(handle, "@%s %s\t%s\tFP:i:%c\tRV:i:%i\tFM:i:%i\n%s\n+\n%s\n", kfp->barcode,
			bufs->FABuffer, bufs->PVBuffer,
			kfp->pass_fail, kfp->n_rc, kfp->length,
			bufs->cons_seq_buffer, kfp->max_phreds);
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
#if DBG
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
	ret.nuc_counts = (uint16_t *)calloc(readlen * 5, sizeof(uint16_t));
	ret.phred_sums = (uint32_t *)calloc(readlen * 4, sizeof(uint32_t));
	ret.max_phreds = (char *)calloc(readlen + 1, sizeof(char)); // Keep track of the maximum phred score observed at position.
	return;
}


static inline KingFisher_t *init_kfp(size_t readlen)
{
	KingFisher_t *ret = (KingFisher_t *)malloc(sizeof(KingFisher_t));
	ret->length = 0; // Check to see if this is necessary after calloc - I'm pretty sure not.
	ret->n_rc = 0;
	ret->readlen = readlen;
	ret->max_phreds = (char *)calloc(readlen + 1, sizeof(char)); // Keep track of the maximum phred score observed at position.
	ret->nuc_counts = (uint16_t *)calloc(readlen * 5, sizeof(uint16_t));
	ret->phred_sums = (uint32_t *)calloc(readlen * 4, sizeof(uint32_t));
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
	char rv;
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

static inline int bc_flip(char *barcode, int blen)
{
	switch(lex_lt(barcode, blen)) {
	case 0:
		return 0;
	case 1:
		return 1;
	default: return rclex_lt(barcode, blen);
	}
}


static inline void pair_flip(mseq_t *m1, mseq_t *m2, tmpvars_t *tmp)
{
	// Test to see if I should flip.
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

static inline void mseq2fq_inline(FILE *handle, mseq_t *mvar, char pass_fail, char *barcode)
{
	fprintf(handle, "@%s ~#!#~|FP=%c|BS=%s|RV=%c\n%s\n+\n%s\n",
			mvar->name, pass_fail, barcode, mvar->rv, mvar->seq, mvar->qual);
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
static inline mseq_t *p7_mseq_rescale_init(kseq_t *seq, char *rescaler, int is_read2)
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
	ret->l = seq->seq.l;
	memset(ret->barcode, 0, MAX_BARCODE_LENGTH + 1);
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
 */
static inline mseq_t *mseq_rescale_init(kseq_t *seq, char *rescaler, tmp_mseq_t *tmp, int is_read2)
{
	mseq_t *ret = p7_mseq_rescale_init(seq, rescaler, is_read2);
	//fprintf(stderr, "Pointer to ret: %p. To tmp: %p. Barcode: %s.\n", ret, tmp, ret->barcode);
	if(!tmp) {
		fprintf(stderr, "Tmpvars not allocated!\n");
		exit(EXIT_FAILURE);
	}
	ret->blen = tmp->blen;
	ret->rv = '0';
	return ret;
}

/*
 * Set is_read2 to 1 for read 2, 0 for read 1.
 */
static inline void update_mseq(mseq_t *mvar, kseq_t *seq, char *rescaler, tmp_mseq_t *tmp, int n_len, int is_read2, int switch_reads)
{
	memcpy(mvar->name, seq->name.s, seq->name.l);
    mvar->name[seq->name.l] = '\0';
	memcpy(mvar->seq, seq->seq.s, seq->seq.l * sizeof(char));
	memset(mvar->seq, 'N', n_len);
	memset(mvar->qual, '#', n_len);
	if(!rescaler) {
		memcpy(mvar->qual + n_len, seq->qual.s + n_len, seq->qual.l * sizeof(char) - n_len);
	}
	else {
		for(uint32_t i = n_len; i < seq->seq.l; i++) {
			// Leave quality scores alone for bases which are N. Otherwise
			mvar->qual[i] = (mvar->seq[i] == 'N') ? '#' : rescale_qscore(is_read2, seq->qual.s[i], i, mvar->seq[i], seq->seq.l, rescaler);
		}
	}
#if DBG
	if(strlen(mvar->qual) != seq->qual.l){
		fprintf(stderr, "Ret qual has the wrong length. (%"PRIu64"). Expected: %"PRIu64". Seq: %s. Kseq: %s.\n", strlen(mvar->qual), seq->qual.l, mvar->qual, seq->qual.s);
		exit(1);
	}
#endif
	mvar->rv = switch_reads ? '1': '0';
}


static inline mseq_t *init_crms_mseq(kseq_t *seq, char *barcode, char *rescaler, tmp_mseq_t *tmp, int is_read2)
{
#if DBG && NDEBEG
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
	mseq_t *ret = mseq_rescale_init(seq, rescaler, tmp, is_read2);
	strcpy(ret->barcode, barcode);
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

/*
 * :param: kseq_t *seq1 - fastq kseq handle
 * :param: kseq_t *seq2 - fastq kseq handle
 * :param: char *barcode - buffer set by function
 * :returns: int - whether or not to switch
 */
static inline int switch_test(kseq_t *seq1, kseq_t *seq2, int offset)
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
	if(lex_strlt(seq1->seq.s + offset, seq2->seq.s + offset)) { // seq1's barcode is lower. No switching.
		memcpy(barcode, seq1->seq.s + offset, blen1_2 * sizeof(char)); // Copying the first half of the barcode
		memcpy(barcode + blen1_2, seq2->seq.s + offset,
				blen1_2 * sizeof(char));
		barcode[blen1_2 * 2] = '\0';
		return 0;
	}
	else {
		memcpy(barcode, seq2->seq.s + offset, blen1_2 * sizeof(char)); // Copying the first half of the barcode
		memcpy(barcode + blen1_2, seq1->seq.s + offset,
				blen1_2 * sizeof(char));
		barcode[blen1_2 * 2] = '\0';
		return 1;
	}
}

static void print_kf(KingFisher_t *kfp)
{
	fprintf(stderr, "Length: %i.\n", kfp->length);
	fprintf(stderr, "FA: ");
	for(int i = 0; i < kfp->readlen; ++i) {
		fprintf(stderr, ",%i\t", (int)kfp->nuc_counts[i]);
	}
	fprintf(stderr, "\n");
	return;
}


static inline void pushback_kseq(KingFisher_t *kfp, kseq_t *seq, int *nuc_indices, int blen)
{
#if !NDEBUG
	fprintf(stderr, "%p: kfp.%p: seq. %p: nuc_indices, %i.\n", kfp, seq, nuc_indices, blen);
	fprintf(stderr, "View: %s, %s.\n", barcode_mem_view(seq), kfp->barcode);
#endif
	if(!kfp->length) {
		memcpy(kfp->barcode, seq->comment.s + 14, blen);
		kfp->barcode[blen] = '\0';
	}
	for(int i = 0; i < kfp->readlen; i++) {
		nuc_to_pos((seq->seq.s[i]), nuc_indices);
		++kfp->nuc_counts[i * 4 + nuc_indices[1]];
#if !NDEBUG
		if(kfp->nuc_counts[i * 4 + nuc_indices[1]] > kfp->length) {
			fprintf(stderr, "Warning: KF counts is too high!\n");
			print_kf(kfp);
			exit(EXIT_FAILURE);
		}
#endif
		kfp->phred_sums[i * 4 + nuc_indices[0]] += seq->qual.s[i] - 33;
		if(seq->qual.s[i] > kfp->max_phreds[i]) {
			kfp->max_phreds[i] = seq->qual.s[i];
		}
	}
	kfp->n_rc += *(barcode_mem_view(seq) + blen + 4) - '0'; // Convert to int
	++kfp->length; // Increment
	return;
}

#endif
