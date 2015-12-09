#ifndef KINGFISHER_H
#define KINGFISHER_H

#include "kseq.h"
#include "stdio.h"
#include "math.h"
#include "charcmp.h"
#include "cstr_util.h"
#include "khash.h"
#include "mem_util.h"
#include "igamc_cephes.h"
#include <zlib.h>
#include <inttypes.h>
#include <assert.h>

#ifndef MAX_BARCODE_LENGTH
#define MAX_BARCODE_LENGTH 30
#endif

#ifndef KSEQ_DEC_GZ
#define KSEQ_DEC_GZ
KSEQ_INIT(gzFile, gzread)
#endif

#define ARGMAX_STR "ACGTN"
#define ARRG_MAX_TO_NUC(argmaxret) ARGMAX_STR[argmaxret]

typedef struct tmpbuffers {
	char name_buffer[120];
	char PVBuffer[1000];
	char FABuffer[1000];
	char cons_seq_buffer[SEQBUF_SIZE];
	uint32_t cons_quals[SEQBUF_SIZE];
	uint16_t agrees[SEQBUF_SIZE];
} tmpbuffers_t;

typedef struct tmpvars {
	char *bs_ptr;
	int blen;
	int readlen;
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
} KingFisher_t;


extern double igamc(double a, double x);
static inline char rescale_qscore(int readnum, int qscore, int cycle, char base, int readlen, char *rescaler);
void destroy_kf(KingFisher_t *kfp);
void stranded_process_write(KingFisher_t *kfpf, KingFisher_t *kfpr, FILE *handle, tmpbuffers_t *bufs);


CONST static inline int ARRG_MAX(KingFisher_t *kfp, int index)
{
	const int i5 = index * 5;
	if(kfp->phred_sums[i5] > kfp->phred_sums[i5 + 1] &&
		kfp->phred_sums[i5] > kfp->phred_sums[i5 + 2] &&
		kfp->phred_sums[i5] > kfp->phred_sums[i5 + 3] &&
		kfp->phred_sums[i5] > kfp->phred_sums[i5 + 4])
		return 0;
	else if(kfp->phred_sums[i5 + 1] > kfp->phred_sums[i5 + 2] &&
			kfp->phred_sums[i5 + 1] > kfp->phred_sums[i5 + 3] &&
			kfp->phred_sums[i5 + 1] > kfp->phred_sums[i5 + 4])
		return 1;
	else if(kfp->phred_sums[i5 + 2] > kfp->phred_sums[i5 + 3] &&
			kfp->phred_sums[i5 + 2] > kfp->phred_sums[i5 + 4])
		return 2;
	else if(kfp->phred_sums[i5 + 3] > kfp->phred_sums[i5 + 4])
		return 3;
	return 4; // 'N'
}


static inline void fill_pv_buffer(int readlen, uint32_t *phred_values, char *buffer)
{
	fill_csv_buffer(readlen, phred_values, buffer, "PV:B:I");
}

static inline void fill_fa_buffer(int readlen, uint16_t *agrees, char *buffer)
{
	char tmpbuf[7];
	strcpy(buffer, "FA:B:I");
	for(int i = 0; i < readlen; ++i) {
		sprintf(tmpbuf, ",%" PRIu16 "", agrees[i]);
		strcat(buffer, tmpbuf);
	}
	return;
}


void dmp_process_write(KingFisher_t *kfp, FILE *handle, tmpbuffers_t *bufs);


CONST static inline char rescale_qscore(int readnum, int qscore, int cycle, char base, int readlen, char *rescaler)
{
	int index = readnum;
	int mult = 2;
	//fprintf(stderr, "index value is now: %i, mult %i.\n", index, mult);
	index += cycle * mult;
	mult *= readlen;
	//fprintf(stderr, "index value is now: %i, mult %i. Qscore: %i, Qscore index%i.\n", index, mult, qscore, qscore - 35);
	index += (qscore - 35) * mult; // Subtract 35 - 33 to get to phred space, 2 to offset by 2.
	mult *= nqscores;
	//fprintf(stderr, "index value is now: %i, mult %i.\n", index, mult);
	index += mult * nuc2num(base);
	//fprintf(stderr, "Index = %i.\n", index);
#if DBG
	if(index >= readlen * 2 * nqscores * 4 || index < 0) {
		fprintf(stderr, "Something's wrong. Index (%i) is too big or negative! Max: %i.\n", index, readlen * 2 * nqscores * 4);
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


KingFisher_t *init_kfp(size_t readlen);


// mseq is a mutable struct holding kseq's information.

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


/*
 * Warning: returns a NULL upon not finding a second pipe symbol.
 * This is *NOT* a properly null-terminated string.
 */
CONST static inline char *barcode_mem_view(kseq_t *seq)
{
	int hits = 0;
	for(int i = 0; i < seq->comment.l; ++i) {
		if(seq->comment.s[i] == '|' || seq->comment.s[i] == '\0') {
			if(!hits) ++hits;
			else
				return (char *)(seq->comment.s + i + 4); // 4 for "|BS="
		}
	}
	return NULL;
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


static inline void tmp_mseq_destroy(tmp_mseq_t mvar)
{
	cond_free(mvar.tmp_seq);
	cond_free(mvar.tmp_qual);
	cond_free(mvar.tmp_barcode);
	mvar.readlen = 0;
	mvar.blen = 0;
}

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
	mseq_t *ret = (mseq_t *)calloc(1, sizeof(mseq_t));
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
	memset(mvar->seq, 'N', n_len), memset(mvar->qual, '#', n_len);
	if(!rescaler)
		memcpy(mvar->qual + n_len, seq->qual.s + n_len, seq->qual.l * sizeof(char) - n_len);
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
}


static inline mseq_t *init_crms_mseq(kseq_t *seq, char *barcode, char *rescaler, tmp_mseq_t *tmp, int is_read2)
{
#if !NDEBEG
	if(!barcode) {
		fprintf(stderr, "Barocde is NULL ABORT>asdfasfjafjhaksdfkjasdfas.\n");
		exit(1);
	}
	else {
		fprintf(stderr, "Barcode: %s.\n", barcode);
	}
	if(rescaler) {
		for(int i = 0; i < seq->seq.l * nqscores * 4 * 2; ++i) {
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


static inline void pb_pos(KingFisher_t *kfp, kseq_t *seq, int i) {
	const uint32_t posdata = nuc2num(seq->seq.s[i]);
	++kfp->nuc_counts[i * 5 + posdata];
	kfp->phred_sums[i * 5 + posdata] += seq->qual.s[i] - 33;
	if(seq->qual.s[i] > kfp->max_phreds[i])
		kfp->max_phreds[i] = seq->qual.s[i];
}

static inline void pushback_kseq(KingFisher_t *kfp, kseq_t *seq, int blen)
{
	if(!(kfp->length++)) { // Increment while checking
		memcpy(kfp->barcode, seq->comment.s + 14, blen);
		kfp->barcode[blen] = '\0';
	}
	for(int i = 0; i < kfp->readlen; ++i) {
		pb_pos(kfp, seq, i);
	}
	return;
}

#endif /*KINGFISHER_H*/
