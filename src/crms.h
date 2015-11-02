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
#include "cstr_utils.h"

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
#define CHECK_CALL(buff, ret) \
	fprintf(stderr, "Now check calling command '%s'.\n", buff); \
	ret = system(buff);\
	if(ret < 0)\
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

static inline int test_homing_seq(kseq_t *seq1, kseq_t *seq2, mssi_settings_t *settings_ptr)
{
	if(!settings_ptr->homing_sequence) {
		return 1;
	}
	else {
		return memcmp(seq1->seq.s + (settings_ptr->blen / 2 + settings_ptr->offset),
					   settings_ptr->homing_sequence,
					   settings_ptr->homing_sequence_length) == 0;
	}
}



static inline char test_hp_inline(char *barcode, int length, int threshold)
{
	int run = 0;
	char last = '\0';
	for(int i = 0; i < length; i++){
		if(barcode[i] == 'N') {
			return '0';
		}
		if(barcode[i] == last) {
			run += 1;
		}
		else {
			run = 0;
			last = barcode[i];
		}
	}
	return (run < threshold) ? '1': '0';
}



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

// Inline function declarations
/*
mseq_t *init_crms_mseq(kseq_t *seq, char *barcode, char *rescaler, tmp_mseq_t *tmp, int n_len, int is_read2);
int crc_flip(mseq_t *mvar, char *barcode, int blen, int readlen);
int ipow(int base, int exp);
void mseq_destroy(mseq_t *mvar);
mseq_t *mseq_rescale_init(kseq_t *seq, char *rescaler, tmp_mseq_t *tmp, int n_len, int is_read2);
int nuc2num(char character);
int nuc_cmp(char forward, char reverse);
char ****parse_rescaler(char *qual_rescale_fname);
char *parse_1d_rescaler(char *qual_rescale_fname);
char rescale_qscore(int readnum, int qscore, int cycle, char base, int readlen, char *rescaler);
//void set_barcode(kseq_t *seq1, kseq_t *seq2, char *barcode, int offset, int blen1_2);
int test_homing_seq(kseq_t *seq1, kseq_t *seq2, mssi_settings_t *settings_ptr);
void tmp_mseq_destroy(tmp_mseq_t mvar);
char nuc_cmpl(char character);
char *make_default_outfname(char *fname, const char *suffix);
char *make_crms_outfname(char *fname);
int get_fileno_limit();
void increase_nofile_limit(int new_limit);
uint64_t get_binnerul(char *barcode, int length);
int get_binner(char *barcode, int length);
uint64_t ulpow(uint64_t base, uint64_t exp);

static splitterhash_params_t *init_splitterhash(mssi_settings_t *settings_ptr, mark_splitter_t *splitter_ptr);
static inline splitterhash_params_t *init_splitterhash_mss(mss_settings_t *settings_ptr, mark_splitter_t *splitter_ptr);
int vl_homing_loc(kseq_t *seq1, kseq_t *seq2, crms_settings_t *settings_ptr);
blens_t *get_blens(char *str2parse);
void free_crms_settings(crms_settings_t settings);
splitterhash_params_t *init_vl_splitterhash(crms_settings_t *settings_ptr, mark_splitter_t *splitter_ptr);
void hash_dmp_core(char *infname, char *outfname);
mseq_t *p7_mseq_rescale_init(kseq_t *seq, char *rescaler, int n_len, int is_read2);
*/

static void splitterhash_destroy(splitterhash_params_t *params);

static inline void FREE_SPLITTER(mark_splitter_t var)
{
	for(int i = 0; i < var.n_handles; i++) {
		//fprintf(stderr, "Now trying to close file #%i with filename %s.\n", i, var.fnames_r1[i]);
		//fprintf(stderr, "Now trying to access FILE * with number %i.\n", i);
		if(var.fnames_r1[i]) {
			free(var.fnames_r1[i]);
			var.fnames_r1[i] = NULL;
		}
		if(var.fnames_r2[i]) {
			free(var.fnames_r2[i]);
			var.fnames_r2[i] = NULL;
		}
	}
	free(var.tmp_out_handles_r1);
	free(var.tmp_out_handles_r2);
	return;
}


static inline void FREE_SPLITTER_PTR(mark_splitter_t *var)
{
	for(int i = 0; i < var->n_handles; i++) {
		//fprintf(stderr, "Now trying to close file #%i with filename %s.\n", i, var->fnames_r1[i]);
		//fprintf(stderr, "Now trying to access FILE * with number %i.\n", i);
		if(var->fnames_r1[i]) {
			free(var->fnames_r1[i]);
			var->fnames_r1[i] = NULL;
		}
		if(var->fnames_r2[i]) {
			free(var->fnames_r2[i]);
			var->fnames_r2[i] = NULL;
		}
	}
	free(var->tmp_out_handles_r1);
	free(var->tmp_out_handles_r2);
	free(var);
	var = NULL;
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


static mark_splitter_t *splitmark_core_rescale(mss_settings_t *settings)
{
	//fprintf(stderr, "[splitmark_rescale_core]\n");
	if(strcmp(settings->input_r1_path, settings->input_r2_path) == 0) {
		fprintf(stderr, "Input read paths are the same {'R1': %s, 'R2': %s}. WTF!\n", settings->input_r1_path, settings->input_r2_path);
		exit(EXIT_FAILURE);
	}
	else {
		fprintf(stderr, "Path to index fq: %s.\n", settings->index_fq_path);
	}
	gzFile fp_read1, fp_read2, fp_index;
	kseq_t *seq1, *seq2, *seq_index;
	mseq_t *rseq1, *rseq2;
	int l1, l2, l_index;
	//fprintf(stderr, "[splitmark_rescale_core]: initializing splitter.\n");
	mark_splitter_t *splitter_ptr = (mark_splitter_t *)malloc(sizeof(mark_splitter_t));
	*splitter_ptr = init_splitter(settings);
	//fprintf(stderr, "[splitmark_rescale_core]: Opening input handles.\n");
	fp_read1 = gzopen(settings->input_r1_path, "r");
	fp_read2 = gzopen(settings->input_r2_path, "r");
	fp_index = gzopen(settings->index_fq_path, "r");
	seq1 = kseq_init(fp_read1);
	seq2 = kseq_init(fp_read2);
	seq_index = kseq_init(fp_index);
	//fprintf(stderr, "[splitmark_rescale_core]: Reading from Read 1. Path: %s.\n", settings->input_r1_path);
	l1 = kseq_read(seq1);
	//fprintf(stderr, "[splitmark_rescale_core]: Reading from Read index. Path: %s.\n", settings->index_fq_path);
	l_index = kseq_read(seq_index);
	//fprintf(stderr, "[splitmark_rescale_core]: Reading from Read 2. Path: %s.\n", settings->input_r2_path);
	l2 = kseq_read(seq2);
	//fprintf(stderr, "[splitmark_rescale_core]: Read from Read 2! Path: %s.\n", settings->input_r2_path);
	int bin = 0;
	int count = 0;
	char pass_fail = '1';
	tmp_mseq_t *tmp = init_tm_ptr(seq1->seq.l, seq_index->seq.l + 2 * settings->salt);
	if(l1 < 0 || l2 < 0 || l_index < 0) {
		fprintf(stderr, "Could not read input fastqs. Abort mission!\n");
		exit(EXIT_FAILURE);
	}
	fprintf(stderr, "Splitter now opening files R1 ('%s'), R2 ('%s'), index ('%s').\n",
			settings->input_r1_path, settings->input_r2_path, settings->index_fq_path);
	rseq1 = p7_mseq_rescale_init(seq1, settings->rescaler, 0); // rseq1 is initialized
	rseq2 = p7_mseq_rescale_init(seq2, settings->rescaler, 1); // rseq2 is initialized
	memcpy(rseq1->barcode, seq1->seq.s + settings->offset, settings->salt); // Copy in the appropriate nucleotides.
	memcpy(rseq1->barcode + settings->salt, seq_index->seq.s, seq_index->seq.l); // Copy in the barcode
	memcpy(rseq1->barcode + settings->salt + seq_index->seq.l, seq2->seq.s + settings->offset, settings->salt);
	rseq1->barcode[settings->salt * 2 + seq_index->seq.l] = '\0';
	update_mseq(rseq1, seq1, settings->rescaler, tmp, 0, 0, 0);
	update_mseq(rseq2, seq2, settings->rescaler, tmp, 0, 1, 0);
	pass_fail = test_hp(rseq1->barcode, settings->hp_threshold);
	bin = get_binner(rseq1->barcode, settings->n_nucs);
	SALTED_MSEQ_2_FQ(splitter_ptr->tmp_out_handles_r1[bin], rseq1, rseq1->barcode, pass_fail);
	SALTED_MSEQ_2_FQ(splitter_ptr->tmp_out_handles_r2[bin], rseq2, rseq1->barcode, pass_fail);
	//fprintf(stderr, "Now beginning to loop through file.\n");
	while ((l1 = kseq_read(seq1)) >= 0 && (l2 = kseq_read(seq2) >= 0)
			&& (l_index = kseq_read(seq_index)) >= 0) {
		if(!(++count % settings->notification_interval)) {
			fprintf(stderr, "Number of records processed: %i.\n", count);
		}
		memcpy(rseq1->barcode, seq1->seq.s + settings->offset, settings->salt); // Copy in the appropriate nucleotides.
		memcpy(rseq1->barcode + settings->salt, seq_index->seq.s, seq_index->seq.l); // Copy in the barcode
		memcpy(rseq1->barcode + settings->salt + seq_index->seq.l, seq2->seq.s + settings->offset, settings->salt);
		update_mseq(rseq1, seq1, settings->rescaler, tmp, 0, 0, 0);
		update_mseq(rseq2, seq2, settings->rescaler, tmp, 0, 1, 0);
		pass_fail = test_hp(rseq1->barcode, settings->hp_threshold);
		bin = get_binner(rseq1->barcode, settings->n_nucs);
		SALTED_MSEQ_2_FQ(splitter_ptr->tmp_out_handles_r1[bin], rseq1, rseq1->barcode, pass_fail);
		SALTED_MSEQ_2_FQ(splitter_ptr->tmp_out_handles_r2[bin], rseq2, rseq1->barcode, pass_fail);
	}
	tm_destroy(tmp);
	mseq_destroy(rseq1);
	mseq_destroy(rseq2);
	kseq_destroy(seq1);
	kseq_destroy(seq2);
	kseq_destroy(seq_index);
	gzclose(fp_read1);
	gzclose(fp_read2);
	gzclose(fp_index);
	for(count = 0; count < settings->n_handles; ++count) {
		fclose(splitter_ptr->tmp_out_handles_r1[count]);
		fclose(splitter_ptr->tmp_out_handles_r2[count]);
		splitter_ptr->tmp_out_handles_r1[count] = NULL;
		splitter_ptr->tmp_out_handles_r2[count] = NULL;
	}
	return splitter_ptr;
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

static inline int nlen_homing_default(kseq_t *seq1, kseq_t *seq2, mssi_settings_t *settings_ptr, int default_len, char pass_fail)
{
	if(settings_ptr->max_blen < 0) {
		pass_fail = (memcmp(seq1->seq.s + (settings_ptr->blen1_2 + settings_ptr->offset),
				   settings_ptr->homing_sequence,
				   settings_ptr->homing_sequence_length) == 0) ? '1': '0';
		return default_len;
	}
	for(int i = settings_ptr->blen1_2 + settings_ptr->offset; i <= settings_ptr->max_blen; ++i) {
		if(memcmp(seq1->seq.s, settings_ptr->homing_sequence, settings_ptr->homing_sequence_length) == 0) {
			pass_fail = '1';
			return i + settings_ptr->homing_sequence_length;
		}
	}
	pass_fail = '0';
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
