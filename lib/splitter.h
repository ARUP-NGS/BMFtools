#ifndef SPLITTER_H
#define SPLITTER_H

#include <stdlib.h>
#include <string.h>
#include "dlib/compiler_util.h"
#include "dlib/mem_util.h"
#include "lib/binner.h"

#ifndef METASYNTACTIC_FNAME_BUFLEN
#define METASYNTACTIC_FNAME_BUFLEN 100
#endif

typedef struct marksplit_settings {
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
	char *index_fq_path; // Make sure this is null if it's inline!
	int max_blen;
	int n_handles; // Number of handles
	int n_nucs; // Number of nucleotides to split by.
	int notification_interval; // How many sets of records do you want to process between progress reports?
	int offset; // Number of bases at the start of the inline barcodes to skip for low quality.
	char *tmp_basename;
	int panthera;
	char *rescaler; // Four-dimensional rescaler array. Size: [readlen, nqscores, 4] (length of reads, number of original quality scores, number of bases)
	char *rescaler_path; // Path to rescaler for
	int run_hash_dmp;
	int salt;
	int threads;
	int is_se;
	int to_stdout;
} marksplit_settings_t;

#ifdef __cplusplus
extern "C" {
#endif
void free_marksplit_settings_ptr(marksplit_settings_t *settings);
void free_marksplit_settings(marksplit_settings_t settings);
#ifdef __cplusplus
}
#endif

typedef struct mark_splitter {
	FILE **tmp_out_handles_r1;
	FILE **tmp_out_handles_r2;
	int n_nucs;
	int n_handles;
	char **fnames_r1;
	char **fnames_r2;
} mark_splitter_t;

#ifdef __cplusplus
extern "C" {
#endif
mark_splitter_t init_splitter(marksplit_settings_t* settings_ptr);
void splitter_destroy(mark_splitter_t *var);
#ifdef __cplusplus
}
#endif

typedef struct splitterhash_params {
	char **infnames_r1;
	char **infnames_r2;
	char **outfnames_r1;
	char **outfnames_r2;
	int n; // Number of infnames and outfnames
	int paired; // 1 if paired, 0 if single-end
} splitterhash_params_t;

#ifdef __cplusplus
extern "C" {
#endif
void splitterhash_destroy(splitterhash_params_t *params);
splitterhash_params_t *init_splitterhash(marksplit_settings_t *settings_ptr, mark_splitter_t *splitter_ptr);
#ifdef __cplusplus
}
#endif

#endif /* SPLITTER_H */
