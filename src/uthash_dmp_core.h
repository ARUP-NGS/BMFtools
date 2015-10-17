#ifndef UTHASH_DMP_CORE_H
#define UTHASH_DMP_CORE_H

#include "crms.h"

void omgz_core(char *infname, char *outfname);


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

#endif
