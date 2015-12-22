#ifndef CRMS_H
#define CRMS_H

#include <getopt.h>
#include <math.h>
#include <omp.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include "binner.h"
#include "cstr_util.h"
#include "kingfisher.h"
#include "io_util.h"
#include "mem_util.h"
#include "nix_util.h"
#include "bmf_hashdmp.h"

#ifndef MAX_HOMING_SEQUENCE
#define MAX_HOMING_SEQUENCE 8
#endif
#ifndef CAT_BUFFER_SIZE
#define CAT_BUFFER_SIZE 500000
#endif
#ifndef MAX_N_BLENS
#define MAX_N_BLENS 10
#endif


#ifndef MAX_BARCODE_LENGTH
#define MAX_BARCODE_LENGTH 30
#endif
#ifndef KSEQ_2_FQ
#define KSEQ_2_FQ(handle, read, index, pass_fail) fprintf(handle, \
		"@%s ~#!#~|FP=%c|BS=%s\n%s\n+\n%s\n",\
	read->name.s, pass_fail + '0', index->seq.s, read->seq.s, read->qual.s)
#endif
#ifndef SALTED_KSEQ_2_FQ
#define SALTED_KSEQ_2_FQ(handle, read, barcode, pass_fail) fprintf(handle, \
		"@%s ~#!#~|FP=%c|BS=%s|RV=0\n%s\n+\n%s\n",\
	read->name.s, pass_fail + '0', barcode, read->seq.s, read->qual.s)
#endif
#ifndef SALTED_MSEQ_2_FQ
#define SALTED_MSEQ_2_FQ(handle, read, barcode, pass_fail) \
	fprintf(handle, \
		"@%s ~#!#~|FP=%c|BS=%s|RV=0\n%s\n+\n%s\n",\
	read->name, pass_fail + '0', barcode, read->seq, read->qual)
#endif


#ifndef FREE_SETTINGS
#define FREE_SETTINGS(settings) free(settings.tmp_basename);\
	free(settings.index_fq_path);\
	free(settings.input_r1_path);\
	free(settings.input_r2_path)
#endif

#ifndef METASYNTACTIC_FNAME_BUFLEN
#define METASYNTACTIC_FNAME_BUFLEN 100
#endif

#ifndef CHECK_CALL
#if !NDEBUG
#define CHECK_CALL(buff) \
	fprintf(stderr, "[D:%s] Now check calling command '%s'.\n", __func__, buff); \
	if(system(buff) < 0) fprintf(stderr, "[D:%s] System call failed. Command: '%s'.\n", __func__, buff)
#else
#define CHECK_CALL(buff) \
	if(system(buff) < 0) fprintf(stderr, "[D:%s] System call failed. Command: '%s'.\n", __func__, buff)
#endif
#endif

extern void hash_dmp_core(char *infname, char *outfname);
extern int isfile(char *);


typedef struct blens {
	int max_blen; // Last value in blens
	int min_blen; // Lowest value in blens
	int blens[MAX_N_BLENS]; // Array holding blens
	int n; // Number of blens to look for
	int current_blen;
} blens_t;

typedef struct crms_settings {
	int hp_threshold; // The minimum length of a homopolymer run to fail a barcode.
	int n_nucs; // Number of nucleotides to split by.
	char *tmp_basename;
	char *input_r1_path;
	char *input_r2_path;
	int n_handles; // Number of handles
	int notification_interval; // How many sets of records do you want to process between progress reports?
	blens_t *blen_data;
	int offset; // Number of bases at the start of the inline barcodes to skip for low quality.
	char ****rescaler; // Three-dimensional rescaler array. Size: [readlen, nqscores, 4] (length of reads, number of original quality scores, number of bases).
	char *rescaler_path; // Path to flat text file for parsing in the rescaler.
	char *ffq_prefix; // Final fastq prefix.
	int threads; // Number of threads to use for parallel dmp.
	char *homing_sequence;
	int homing_sequence_length;
	int run_hash_dmp;
} crms_settings_t;




//char *trim_ext(char *fname);


int test_homing_seq(kseq_t *seq1, kseq_t *seq2, marksplit_settings_t *settings_ptr);
char test_hp_inline(char *barcode, int length, int threshold);
void clean_homing_sequence(char *);



CONST static inline int test_hp(char *barcode, int threshold)
{
	int run = 0; char last = '\0';
	while(*barcode) {
		if(*barcode == 'N') return 0;
		if(*barcode == last) {
			if(++run == threshold)
				return 0;
		} else {
			last = *barcode; run = 0;
		}
		barcode++;
	}
	return 1;
}


/*
 * :param: settings [crms_settings_t, marksplit_settings_t] Settings struct in which to free the rescaler.
 * :return: void
 * This function supersedes free_rescaler_array by being type-generic.
 */

#define cfree_rescaler(settings) \
	do {\
		if(settings.rescaler) {\
			int readlen##_settings = count_lines(settings.rescaler_path);\
			for(int i = 0; i < 2; ++i) {\
				for(int j = 0; j < readlen##_settings; ++j) {\
					for(int k = 0; k < nqscores; ++k) {\
						cond_free(settings.rescaler[i][j][k]);\
					}\
					cond_free(settings.rescaler[i][j]);\
				}\
				cond_free(settings.rescaler[i]);\
			}\
			cond_free(settings.rescaler);\
		}\
	} while(0)


extern splitterhash_params_t *init_splitterhash(marksplit_settings_t *settings_ptr, mark_splitter_t *splitter_ptr);
extern void splitterhash_destroy(splitterhash_params_t *params);

/*
 * Returns a null-terminated string with the default outfname.
 * Warning: Must be freed!
 */
static inline char *make_default_outfname(char *fname, const char *suffix)
{
	char buf[200];
	char *prefix = trim_ext(fname);
	strcpy(buf, prefix);
	strcat(buf, suffix);
	char *ret = strdup(buf);
	free(prefix);
	return ret;
}

static inline char *make_crms_outfname(char *fname)
{
	return make_default_outfname(fname, ".crms.split");
}


CONST static inline int infer_barcode_length(char *bs_ptr)
{
	int ret = 0;
	for (;;++ret)
		if(bs_ptr[ret] == '\0' || bs_ptr[ret] == '|')
			return ret;
	return -1; // This never happens.
}


#ifndef kroundup32
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif


CONST static inline int nlen_homing_seq(kseq_t *seq1, kseq_t *seq2, marksplit_settings_t *settings_ptr)
{
	if(settings_ptr->max_blen < 0) {
		return (memcmp(seq1->seq.s + (settings_ptr->blen1_2 + settings_ptr->offset),
					   settings_ptr->homing_sequence,
					   settings_ptr->homing_sequence_length) == 0) ? settings_ptr->blen1_2 + settings_ptr->offset + settings_ptr->homing_sequence_length: -1;
	}
	for(int i = settings_ptr->blen1_2 + settings_ptr->offset; i <= settings_ptr->max_blen; ++i) {
		if(memcmp(seq1->seq.s, settings_ptr->homing_sequence, settings_ptr->homing_sequence_length) == 0) {
			return i + settings_ptr->homing_sequence_length;
		}
	}
	return -1;
}

CONST static inline int nlen_homing_default(kseq_t *seq1, kseq_t *seq2, marksplit_settings_t *settings_ptr, int default_len, int *pass_fail)
{
	if(settings_ptr->max_blen < 0) {
		*pass_fail = (memcmp(seq1->seq.s + (settings_ptr->blen1_2 + settings_ptr->offset),
				   settings_ptr->homing_sequence,
				   settings_ptr->homing_sequence_length) == 0) ? 1: 0;
		return default_len;
	}
	for(int i = settings_ptr->blen1_2 + settings_ptr->offset; i <= settings_ptr->max_blen; ++i) {
		if(memcmp(seq1->seq.s, settings_ptr->homing_sequence, settings_ptr->homing_sequence_length) == 0) {
			*pass_fail = 1;
			return i + settings_ptr->homing_sequence_length;
		}
	}
	*pass_fail = 0;
	return default_len;
}


#ifndef METASYNTACTIC_FNAME_BUFLEN
#define METASYNTACTIC_FNAME_BUFLEN 100
#endif


int ipow(int base, int exp);


#ifndef KSEQ_2_FQ
#define KSEQ_2_FQ(handle, read, index, pass_fail) fprintf(handle, \
		"@%s ~#!#~|FP=%c|BS=%s\n%s\n+\n%s\n",\
	read->name.s, pass_fail, index->seq.s, read->seq.s, read->qual.s)
#endif


#ifndef FREE_SETTINGS
#define FREE_SETTINGS(settings) cond_free(settings.tmp_basename);\
	cond_free(settings.index_fq_path);\
	cond_free(settings.input_r1_path);\
	cond_free(settings.input_r2_path);\
	cond_free(settings.ffq_prefix)
#endif

#endif
