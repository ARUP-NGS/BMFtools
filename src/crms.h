#ifndef CRMS_H
#define CRMS_H

#include <getopt.h>
#include <math.h>
#include <omp.h>
#include <regex.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/resource.h>
#include <zlib.h>
#include "kingfisher.h"
#include "array_parser.h"
#include "nix_resource.h"
#include "mem_util.h"
#include "binner.h"
#include "cstr_util.h"

#ifndef MAX_HOMING_SEQUENCE
#define MAX_HOMING_SEQUENCE 8
#endif
#ifndef CAT_BUFFER_SIZE
#define CAT_BUFFER_SIZE 500000
#endif
#ifndef METASYNTACTIC_FNAME_BUFLEN
#define METASYNTACTIC_FNAME_BUFLEN 100
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
	read->name.s, pass_fail, index->seq.s, read->seq.s, read->qual.s)
#endif
#ifndef SALTED_KSEQ_2_Q
#define SALTED_KSEQ_2_FQ(handle, read, barcode, pass_fail) fprintf(handle, \
		"@%s ~#!#~|FP=%c|BS=%s|RV=0\n%s\n+\n%s\n",\
	read->name.s, pass_fail, barcode, read->seq.s, read->qual.s)
#endif
#ifndef SALTED_MSEQ_2_FQ
#define SALTED_MSEQ_2_FQ(handle, read, barcode, pass_fail) \
	fprintf(handle, \
		"@%s ~#!#~|FP=%c|BS=%s|RV=0\n%s\n+\n%s\n",\
	read->name, pass_fail, barcode, read->seq, read->qual)
#endif


#ifndef FREE_SETTINGS
#define FREE_SETTINGS(settings) free(settings.output_basename);\
	free(settings.index_fq_path);\
	free(settings.input_r1_path);\
	free(settings.input_r2_path)
#endif

#ifndef METASYNTACTIC_FNAME_BUFLEN
#define METASYNTACTIC_FNAME_BUFLEN 100
#endif

#ifndef CHECK_CALL
#define CHECK_CALL(buff) \
	fprintf(stderr, "Now check calling command '%s'.\n", buff); \
	if(system(buff) < 0)\
		fprintf(stderr, "System call failed. Command: '%s'.\n", buff)
#endif

extern void khash_dmp_core(char *infname, char *outfname);


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
	char *output_basename;
	char *input_r1_path;
	char *input_r2_path;
	int n_handles; // Number of handles
	int notification_interval; // How many sets of records do you want to process between progress reports?
	blens_t *blen_data;
	int offset; // Number of bases at the start of the inline barcodes to skip for low quality.
	char ****rescaler; // Three-dimensional rescaler array. Size: [readlen, 39, 4] (length of reads, number of original quality scores, number of bases).
	char *rescaler_path; // Path to flat text file for parsing in the rescaler.
	char *ffq_prefix; // Final fastq prefix.
	int threads; // Number of threads to use for parallel dmp.
	char *homing_sequence;
	int homing_sequence_length;
	int run_hash_dmp;
} crms_settings_t;



typedef struct mssi_settings {
	int hp_threshold; // The minimum length of a homopolymer run to fail a barcode.
	int n_nucs; // Number of nucleotides to split by.
	char *output_basename;
	char *input_r1_path;
	char *input_r2_path;
	char *homing_sequence; // Homing sequence...
	int homing_sequence_length; // Length of homing sequence, should it be used.
	int n_handles; // Number of handles
	int notification_interval; // How many sets of records do you want to process between progress reports?
	int offset; // Number of bases at the start of the inline barcodes to skip for low quality.
	char *rescaler_path; // Path to rescaler for
	int threads;
	int run_hash_dmp;
	char *ffq_prefix; // Final fastq prefix
	char *rescaler; // Four-dimensional rescaler array. Size: [readlen, 39, 4] (length of reads, number of original quality scores, number of bases)
	int blen;
	int blen1_2;
	int max_blen;
	int gzip_output;
	int panthera;
	int gzip_compression;
	int cleanup; // Set to false to leave temporary files
	int annealed; // Set to true to avoid reversing sequences TODO: Actually implement this.
} mssi_settings_t;


typedef struct mss_settings {
	int cleanup;
	int hp_threshold;
	int n_nucs;
	int n_handles;
	char *output_basename;
	int notification_interval;
	char *index_fq_path;
	char *input_r1_path;
	char *input_r2_path;
	char *ffq_prefix;
	int run_hash_dmp;
	int gzip_output;
	int gzip_compression;
	int panthera; // One "big cat" or many small cats?
	int salt; // Number of bases from each of read 1 and read 2 to use to salt
	int offset; // The number of bases at the start of reads 1 and 2 to skip when salting
	int threads; // Number of threads to use for parallel dmp
	char *rescaler_path;
	char *rescaler;
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


//char *trim_ext(char *fname);


int test_homing_seq(kseq_t *seq1, kseq_t *seq2, mssi_settings_t *settings_ptr);
char test_hp_inline(char *barcode, int length, int threshold);
extern void splitterhash_destroy(splitterhash_params_t *params);



static inline char test_hp(char *seq, int threshold)
{
	int run = 0;
	char last = '\0';
	for(int i = 0; seq[i]; ++i){
		if(seq[i] == 'N') {
			return '0';
		}
		if(seq[i] == last) {
			++run;
		}
		else {
			run = 0;
			last = seq[i];
		}
	}
	return (run < threshold) ? '1': '0';
}


/*
 * :param: settings [crms_settings_t, mssi_settings_t] Settings struct in which to free the rescaler.
 * :return: void
 * This function supersedes free_rescaler_array by being type-generic.
 */

#define cfree_rescaler(settings) \
	do {\
		if(settings.rescaler) {\
			int readlen##_settings = count_lines(settings.rescaler_path);\
			for(int i = 0; i < 2; ++i) {\
				for(int j = 0; j < readlen##_settings; ++j) {\
					for(int k = 0; k < 39; ++k) {\
						cond_free(settings.rescaler[i][j][k]);\
					}\
					cond_free(settings.rescaler[i][j]);\
				}\
				cond_free(settings.rescaler[i]);\
			}\
			cond_free(settings.rescaler);\
		}\
	} while(0)


static inline splitterhash_params_t *init_vl_splitterhash(crms_settings_t *settings_ptr, mark_splitter_t *splitter_ptr)
{
#if DBG
	fprintf(stderr, "Initializing splitterhash. Output basename: %s.\n", settings_ptr->output_basename);
#endif
	if(!settings_ptr) {
		fprintf(stderr, "Settings pointer null. Abort!\n");
		exit(EXIT_FAILURE);
	}
	if(!settings_ptr->output_basename) {
		fprintf(stderr, "Output basename not set. Abort!\n");
		exit(EXIT_FAILURE);
	}
	if(!splitter_ptr) {
		fprintf(stderr, "Splitter pointer null. Abort!\n");
		exit(EXIT_FAILURE);
	}
	char tmp_buffer [METASYNTACTIC_FNAME_BUFLEN];
	splitterhash_params_t *ret = (splitterhash_params_t *)malloc(sizeof(splitterhash_params_t));
	fprintf(stderr, "Alloc'd ret.\n");
	ret->n = splitter_ptr->n_handles;
	ret->outfnames_r1 = (char **)malloc(ret->n * sizeof(char *));
	ret->outfnames_r2 = (char **)malloc(ret->n * sizeof(char *));
	ret->infnames_r1 = (char **)malloc(ret->n * sizeof(char *));
	ret->infnames_r2 = (char **)malloc(ret->n * sizeof(char *));
	for(int i = 0; i < splitter_ptr->n_handles; ++i) {
		ret->infnames_r1[i] = splitter_ptr->fnames_r1[i];
		ret->infnames_r2[i] = splitter_ptr->fnames_r2[i]; // Does not allocate memory.  This is freed by mark_splitter_t!
		sprintf(tmp_buffer, "%s.%i.R1.dmp.fastq", settings_ptr->output_basename, i);
		ret->outfnames_r1[i] = strdup(tmp_buffer);
		sprintf(tmp_buffer, "%s.%i.R2.dmp.fastq", settings_ptr->output_basename, i);
		ret->outfnames_r2[i] = strdup(tmp_buffer);
	}
	return ret;
}


static const char *crms_suffix = ".crms.split";

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
	return make_default_outfname(fname, crms_suffix);
}

#define free_rescaler_array(settings) \
		int readlen = count_lines(settings.rescaler_path);\
		for(int i##settings = 0; i##settings < 2; ++i##settings) {\
			for(int j##settings = 0; j##settings < readlen; ++j##settings) {\
				for(int k##settings = 0; k##settings < 39; ++k##settings) {\
					cond_free(settings.rescaler[i##settings][j##settings][k##settings]);\
				}\
				cond_free(settings.rescaler[i##settings][j##settings]);\
			}\
			cond_free(settings.rescaler[i##settings]);\
		}\
		cond_free(settings.rescaler)\

static inline void FREE_SPLITTER_PTR(mark_splitter_t *var)
{
	for(int i = 0; i < var->n_handles; i++) {
		cond_free(var->fnames_r1[i]); cond_free(var->fnames_r2[i]);
	}
	free(var->tmp_out_handles_r1);
	free(var->tmp_out_handles_r2);
	free(var), var = NULL;
	return;
}

static inline mark_splitter_t init_splitter(mss_settings_t* settings_ptr)
{
	mark_splitter_t ret = {
		.n_handles = ipow(4, settings_ptr->n_nucs),
		.n_nucs = settings_ptr->n_nucs,
		.fnames_r1 = NULL,
		.fnames_r2 = NULL,
		.tmp_out_handles_r1 = NULL,
		.tmp_out_handles_r2 = NULL
	};
	// Avoid passing more variables.
	ret.tmp_out_handles_r1 = (FILE **)malloc(settings_ptr->n_handles * sizeof(FILE *));
	ret.tmp_out_handles_r2 = (FILE **)malloc(settings_ptr->n_handles * sizeof(FILE *));
	ret.fnames_r1 = (char **)malloc(ret.n_handles * sizeof(char *));
	ret.fnames_r2 = (char **)malloc(ret.n_handles * sizeof(char *));
	char tmp_buffer [METASYNTACTIC_FNAME_BUFLEN];
	for (int i = 0; i < ret.n_handles; i++) {
		sprintf(tmp_buffer, "%s.tmp.%i.R1.fastq", settings_ptr->output_basename, i);
		ret.fnames_r1[i] = strdup(tmp_buffer);
		sprintf(tmp_buffer, "%s.tmp.%i.R2.fastq", settings_ptr->output_basename, i);
		ret.fnames_r2[i] = strdup(tmp_buffer);
		ret.tmp_out_handles_r1[i] = fopen(ret.fnames_r1[i], "w");
		ret.tmp_out_handles_r2[i] = fopen(ret.fnames_r2[i], "w");
	}
	return ret;
}


static inline mark_splitter_t init_splitter_inline(mssi_settings_t* settings_ptr)
{
	mark_splitter_t ret = {
		.n_handles = ipow(4, settings_ptr->n_nucs),
		.n_nucs = settings_ptr->n_nucs,
		.fnames_r1 = NULL,
		.fnames_r2 = NULL,
		.tmp_out_handles_r1 = NULL,
		.tmp_out_handles_r2 = NULL
	};
	// Avoid passing more variables.
	ret.tmp_out_handles_r1 = (FILE **)malloc(settings_ptr->n_handles * sizeof(FILE *));
	ret.tmp_out_handles_r2 = (FILE **)malloc(settings_ptr->n_handles * sizeof(FILE *));
	ret.fnames_r1 = (char **)malloc(ret.n_handles * sizeof(char *));
	ret.fnames_r2 = (char **)malloc(ret.n_handles * sizeof(char *));
	char tmp_buffer [METASYNTACTIC_FNAME_BUFLEN];
	for (int i = 0; i < ret.n_handles; i++) {
		sprintf(tmp_buffer, "%s.tmp.%i.R1.fastq", settings_ptr->output_basename, i);
		ret.fnames_r1[i] = strdup(tmp_buffer);
		sprintf(tmp_buffer, "%s.tmp.%i.R2.fastq", settings_ptr->output_basename, i);
		ret.fnames_r2[i] = strdup(tmp_buffer);
		ret.tmp_out_handles_r1[i] = fopen(ret.fnames_r1[i], "w");
		ret.tmp_out_handles_r2[i] = fopen(ret.fnames_r2[i], "w");
	}
	return ret;
}


static inline int infer_barcode_length(char *bs_ptr)
{
	int ret = 0;
	for (;;++ret) {
		if(bs_ptr[ret] == '\0' || bs_ptr[ret] == '|') {
			return ret;
		}
	}
}


static inline void free_mssi_settings(mssi_settings_t settings)
{
	free(settings.output_basename);
	free(settings.input_r1_path);
	free(settings.input_r2_path);
	if(settings.rescaler_path) free(settings.rescaler_path);
	return;
}

#ifndef FREE_MSSI_SETTINGS_PTR
#define FREE_MSSI_SETTINGS_PTR(settings) cond_free(settings->output_basename);\
	cond_free(settings->input_r1_path);\
	cond_free(settings->input_r2_path);\
	cond_free(settings->rescaler_path)
#endif

#ifndef FREE_MSSI_SETTINGS
#define FREE_MSSI_SETTINGS(settings) cond_free(settings.output_basename);\
	cond_free(settings.input_r1_path);\
	cond_free(settings.input_r2_path);\
	cond_free(settings.rescaler_path);
#endif


#ifndef kroundup32
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif


// Calls incomplete gamma complement from CEPHES.

static inline int nlen_homing_seq(kseq_t *seq1, kseq_t *seq2, mssi_settings_t *settings_ptr)
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

static inline int nlen_homing_default(kseq_t *seq1, kseq_t *seq2, mssi_settings_t *settings_ptr, int default_len, char *pass_fail)
{
	if(settings_ptr->max_blen < 0) {
		*pass_fail = (memcmp(seq1->seq.s + (settings_ptr->blen1_2 + settings_ptr->offset),
				   settings_ptr->homing_sequence,
				   settings_ptr->homing_sequence_length) == 0) ? '1': '0';
		return default_len;
	}
	for(int i = settings_ptr->blen1_2 + settings_ptr->offset; i <= settings_ptr->max_blen; ++i) {
		if(memcmp(seq1->seq.s, settings_ptr->homing_sequence, settings_ptr->homing_sequence_length) == 0) {
			*pass_fail = '1';
			return i + settings_ptr->homing_sequence_length;
		}
	}
	*pass_fail = '0';
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
#define FREE_SETTINGS(settings) cond_free(settings.output_basename);\
	cond_free(settings.index_fq_path);\
	cond_free(settings.input_r1_path);\
	cond_free(settings.input_r2_path);\
	cond_free(settings.ffq_prefix)
#endif

#endif
