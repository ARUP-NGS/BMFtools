#ifndef SPLITTER_H
#define SPLITTER_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "compiler_util.h"
#include "mem_util.h"
#include "binner.h"

#ifndef METASYNTACTIC_FNAME_BUFLEN
#define METASYNTACTIC_FNAME_BUFLEN 100
#endif

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

void free_mssi_settings_ptr(mssi_settings_t *settings);
void free_mssi_settings(mssi_settings_t settings);

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

mark_splitter_t init_splitter_inline(mssi_settings_t* settings_ptr);
mark_splitter_t init_splitter(mss_settings_t* settings_ptr);
void splitter_destroy(mark_splitter_t *var);

typedef struct splitterhash_params {
	char **infnames_r1;
	char **infnames_r2;
	char **outfnames_r1;
	char **outfnames_r2;
	int n; // Number of infnames and outfnames
	int paired; // 1 if paired, 0 if single-end
} splitterhash_params_t;

void splitterhash_destroy(splitterhash_params_t *params);
splitterhash_params_t *init_splitterhash(mssi_settings_t *settings_ptr, mark_splitter_t *splitter_ptr);
splitterhash_params_t *init_splitterhash_mss(mss_settings_t *settings_ptr, mark_splitter_t *splitter_ptr);

#endif /* SPLITTER_H */
