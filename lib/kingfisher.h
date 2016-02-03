#ifndef KINGFISHER_H
#define KINGFISHER_H

#include <stdio.h>
#include <math.h>
#include <zlib.h>
#include <inttypes.h>
#include "htslib/khash.h"
#include "htslib/kseq.h"
#include "include/igamc_cephes.h"
#include "lib/mseq.h"
#include "lib/rescaler.h"
#include "lib/splitter.h"
#include "dlib/char_util.h"
#include "dlib/cstr_util.h"
#include "dlib/mem_util.h"
#include "dlib/compiler_util.h"

#define MAX_PV 3117 // Maximum seen with doubles
#define MIN_FRAC_AGREED 0.5 // Minimum fraction of bases agreed in a family to not "N" the base.
#define HASH_DMP_OFFSET 14
#define FP_OFFSET 9

#ifndef KSEQ_DEC_GZ
#define KSEQ_DEC_GZ
KSEQ_INIT(gzFile, gzread)
#endif


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
	char *max_phreds; // Maximum phred score observed at position. Use this as the final sequence for the quality to maintain compatibility with GATK and other tools.
	int length; // Number of reads in family
	int readlen; // Length of reads
	char barcode[MAX_BARCODE_LENGTH + 1];
	char pass_fail;
} KingFisher_t;



extern double igamc(double a, double x);

#ifdef __cpluslus
extern "C" {
#endif


static inline void pushback_kseq(KingFisher_t *kfp, kseq_t *seq, int blen);
static inline void pb_pos(KingFisher_t *kfp, kseq_t *seq, int i);
static inline char rescale_qscore(int readnum, char qscore, int cycle, char base, int readlen, char *rescaler);
void stranded_process_write(KingFisher_t *kfpf, KingFisher_t *kfpr, FILE *handle, tmpbuffers_t *bufs);
void dmp_process_write(KingFisher_t *kfp, FILE *handle, tmpbuffers_t *bufs, int is_rev);
CONST static inline int kfp_argmax(KingFisher_t *kfp, int index);
CONST static inline int arr_max_u32(uint32_t *arr, int index);


/*
 * @func fill_fa
 * Calls append_csv_buffer for 32-bit PV array tags.
 * :param: readlen [int] Length of read
 * :param: arr [uint16_t *] Array of values to put into the buffer.
 * :param: buffer [char *] Buffer for the values.
 */
static inline void fill_fa(int readlen, uint16_t *agrees, char *buffer)
{
	char tmpbuf[7];
	memcpy(buffer, "FA:B:I", 7); // "Copy FA:B:I:\0" over
	for(int i = 0; i < readlen; ++i) {
		sprintf(tmpbuf, ",%" PRIu16 "", agrees[i]);
		strcat(buffer, tmpbuf);
	}
}

static inline void pb_pos(KingFisher_t *kfp, kseq_t *seq, int i) {
	const uint32_t posdata = nuc2num(seq->seq.s[i]) + i * 5;
	++kfp->nuc_counts[posdata];
	kfp->phred_sums[posdata] += seq->qual.s[i] - 33;
	if(seq->qual.s[i] > kfp->max_phreds[posdata]) kfp->max_phreds[posdata] = seq->qual.s[i];
}


static inline void pushback_kseq(KingFisher_t *kfp, kseq_t *seq, int blen)
{
	if(!kfp->length++) { // Increment while checking
		kfp->pass_fail = seq->comment.s[FP_OFFSET];
		memcpy(kfp->barcode, seq->comment.s + HASH_DMP_OFFSET, blen);
		kfp->barcode[blen] = '\0';
	}
	for(int i = 0; i < kfp->readlen; ++i) pb_pos(kfp, seq, i);
}


/*
 * @func arr_max_u32
 * :param: arr [uint32_t *] 2-d array of values. 5 * index + basecall is the index to use.
 * :param: index [int] Base in read to find the maximum value for.
 * :returns: [int] the nucleotide number for the maximum value at this index in the read.
 */
CONST static inline int arr_max_u32(uint32_t *arr, int index)
{
	const uint32_t i5 = index * 5;
	if(arr[i5] > arr[i5 + 1] &&
		arr[i5] > arr[i5 + 2] &&
		arr[i5] > arr[i5 + 3] &&
		arr[i5] > arr[i5 + 4])
		return 0;
	else if(arr[i5 + 1] > arr[i5 + 2] &&
			arr[i5 + 1] > arr[i5 + 3] &&
			arr[i5 + 1] > arr[i5 + 4])
		return 1;
	else if(arr[i5 + 2] > arr[i5 + 3] &&
			arr[i5 + 2] > arr[i5 + 4])
		return 2;
	else if(arr[i5 + 3] > arr[i5 + 4])
		return 3;
	return 4; // 'N'
}


CONST static inline int kfp_argmax(KingFisher_t *kfp, int index)
{
	return arr_max_u32(kfp->phred_sums, index);
}


#ifdef __cpluslus
}
#endif

#endif /*KINGFISHER_H*/
