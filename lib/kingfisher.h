#ifndef KINGFISHER_H
#define KINGFISHER_H

#include "kseq.h"
#include "stdio.h"
#include "math.h"
#include "char_util.h"
#include "cstr_util.h"
#include "khash.h"
#include "mem_util.h"
#include "igamc_cephes.h"
#include "mseq.h"
#include "rescaler.h"
#include <zlib.h>
#include <inttypes.h>
#include <assert.h>

#define MAX_PV 3117 // Maximum seen with doubles
#define MIN_FRAC_AGREED 0.5 // Minimum fraction of bases agreed in a family to not "N" the base.

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


typedef struct mssi_settings {
	int annealed; // Set to true to avoid reversing sequences TODO: Actually implement this.
	int blen;
	int blen1_2;
	int cleanup; // Set to false to leave temporary files
	char *ffq_prefix; // Final fastq prefix
	int gzip_compression;
	int gzip_output;
	char *homing_sequence; // Homing sequence...
	int homing_sequence_length; // Length of homing sequence, should it be used.
	int hp_threshold; // The minimum length of a homopolymer run to fail a barcode.
	char *input_r1_path;
	char *input_r2_path;
	int max_blen;
	int n_handles; // Number of handles
	int n_nucs; // Number of nucleotides to split by.
	int notification_interval; // How many sets of records do you want to process between progress reports?
	int offset; // Number of bases at the start of the inline barcodes to skip for low quality.
	char *output_basename;
	int panthera;
	char *rescaler; // Four-dimensional rescaler array. Size: [readlen, nqscores, 4] (length of reads, number of original quality scores, number of bases)
	char *rescaler_path; // Path to rescaler for
	int run_hash_dmp;
	int threads;
} mssi_settings_t;

typedef struct mss_settings {
	int cleanup;
	char *ffq_prefix;
	int gzip_compression;
	int gzip_output;
	int hp_threshold;
	char *index_fq_path;
	char *input_r1_path;
	char *input_r2_path;
	int n_handles;
	int n_nucs;
	int notification_interval;
	int offset; // The number of bases at the start of reads 1 and 2 to skip when salting
	char *output_basename;
	int panthera; // One "big cat" or many small cats?
	char *rescaler;
	char *rescaler_path;
	int run_hash_dmp;
	int salt; // Number of bases from each of read 1 and read 2 to use to salt
	int threads; // Number of threads to use for parallel dmp
} mss_settings_t;


typedef struct mark_splitter {
	FILE **tmp_out_handles_r1;
	FILE **tmp_out_handles_r2;
	int n_nucs;
	int n_handles;
	char **fnames_r1;
	char **fnames_r2;
} mark_splitter_t;


typedef struct splitterhash_params {
	char **infnames_r1;
	char **infnames_r2;
	char **outfnames_r1;
	char **outfnames_r2;
	int n; // Number of infnames and outfnames
	int paired; // 1 if paired, 0 if single-end
} splitterhash_params_t;


typedef struct KingFisher {
	uint16_t *nuc_counts; // Count of nucleotides of this form
	uint32_t *phred_sums; // Sums of -10log10(p-value)
	int length; // Number of reads in family
	int readlen; // Length of reads
	char *max_phreds; // Maximum phred score observed at position. Use this as the final sequence for the quality to maintain compatibility with GATK and other tools.
	char barcode[MAX_BARCODE_LENGTH + 1];
	char pass_fail;
} KingFisher_t;

KingFisher_t *init_kfp(size_t readlen);
void destroy_kf(KingFisher_t *kfp);
static inline void pushback_kseq(KingFisher_t *kfp, kseq_t *seq, int blen);
static inline void pb_pos(KingFisher_t *kfp, kseq_t *seq, int i);





extern double igamc(double a, double x);


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
	if(seq->qual.s[i] > kfp->max_phreds[posdata])
		kfp->max_phreds[posdata] = seq->qual.s[i];
}

static inline void pushback_kseq(KingFisher_t *kfp, kseq_t *seq, int blen)
{
	if(!kfp->length++) { // Increment while checking
		memcpy(kfp->barcode, seq->comment.s + 14, blen);
		kfp->barcode[blen] = '\0';
	}
	for(int i = 0; i < kfp->readlen; ++i)
		pb_pos(kfp, seq, i);
	return;
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

#endif /*KINGFISHER_H*/
