#pragma once

#include "include/kseq_dec.h"

typedef struct tmpvars {
	char *bs_ptr;
	int blen;
	int readlen;
	int nuc_indices[2];
	char key[MAX_BARCODE_LENGTH + 1];
	int l; // For holding ret value for seq.
	tmpbuffers_t *buffers;
} tmpvars_t;

void omgz_core(char *infname, char *outfname);


typedef struct HashKing {
	UT_hash_handle hh;
	char id[MAX_BARCODE_LENGTH + 1];
	KingFisher_t *value;
} HashKing_t;

typedef struct splitterhash_params {
	char **infnames_r1;
	char **infnames_r2;
	char **outfnames_r1;
	char **outfnames_r2;
	int n; // Number of infnames and outfnames
	int paired; // 1 if paired, 0 if single-end
} splitterhash_params_t;


/*
 * :param: seq - [arg/kseq_t *] a filled-in kseq object.
 * :param: buf - a pre-allocated buffer or malloc'd char_ptr with enough space for the barcode and the null terminus.
 * :returns:
 */
inline void cp_bs2buf(kseq_t *seq, char *buf)
{
	char *view = barcode_mem_view(seq);
	int blen = 0;
	while(view[blen] != '\0' && view[blen] != '|') {
		++blen;
	}
	memcpy(buf, view, blen);
	buf[blen] = '\0';
	return;
}



inline splitterhash_params_t *init_splitterhash_mss(mss_settings_t *settings_ptr, mark_splitter_t *splitter_ptr)
{
#if !NDEBUG
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


inline splitterhash_params_t *init_splitterhash(mssi_settings_t *settings_ptr, mark_splitter_t *splitter_ptr)
{
#if !NDEBUG
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

inline void splitterhash_destroy(splitterhash_params_t *params)
{
	for(int i = 0; i < params->n; ++i) {
		free(params->outfnames_r1[i]);
		free(params->outfnames_r2[i]);
	}
	free(params->outfnames_r1);
	free(params->outfnames_r2);
	free(params);
	params = NULL;
	return;
}


extern double igamc(double a, double x);
extern uint64_t get_binnerul(char *barcode, int length); // From binner.h
extern int64_t get_binnerl(char *barcode, int length); // From binner.h
extern char *barcode_mem_view(kseq_t *seq); // from dmp.h
int64_t lpow(int64_t base, int64_t exp);
char ARRG_MAX_TO_NUC(int argmaxret);
double igamc_pvalues(int num_pvalues, double x);
int ARRG_MAX(KingFisher_t *kfp, int index);
int infer_barcode_length(char *bs_ptr);
int pvalue_to_phred(double pvalue);
void destroy_kf(KingFisher_t *kfp);
void dmp_process_write(KingFisher_t *kfp, FILE *handle, int blen, tmpbuffers_t *tmp);
void fill_csv_buffer(int readlen, int *arr, char *buffer, char *prefix, char typecode);
void fill_fa_buffer(KingFisher_t *kfp, int *agrees, char *buffer);
void fill_pv_buffer(KingFisher_t *kfp, int *agrees, char *buffer);
void pushback_kseq(KingFisher_t *kfp, kseq_t *seq, int *nuc_indices, int blen);
KingFisher_t init_kf(int readlen);
void nuc_to_pos(char character, int *nuc_indices);
char test_hp(char *seq, int threshold);
int get_binner(char *barcode, int length);
void hash_dmp_core(FILE *handle, HashKing_t *hash, tmpvars_t tmp, kseq_t *seq);
uint64_t ulpow(uint64_t base, uint64_t exp);
void cp_view2buf(char *view, char *buf);
void omgz_core(char *infname, char *outfname);
void p7_mseq_rescale_init(kseq_t *seq, mseq_t *ret, char *rescaler, int n_len, int is_read2);
char rescale_qscore(int readnum, int qscore, int cycle, char base, int readlen, char *rescaler);
