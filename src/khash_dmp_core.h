#ifndef KHASH_DMP_CORE_H
#define KHASH_DMP_CORE_H
#include "khash.h"
#include "crms.h"

KHASH_MAP_INIT_STR(dmp, KingFisher_t *)
void khash_dmp_core(char *infname, char *outfname);
int khash_dmp_main(int argc, char *argv[]);

static inline void cp_view2buf(char *view, char *buf)
{
	int blen = 0;
	while(view[blen] != '\0' && view[blen] != '|') {
		++blen;
	}
	memcpy(buf, view, blen);
	buf[blen] = '\0';
	return;
}


static inline tmpvars_t *init_tmpvars_p(char *bs_ptr, int blen, int readlen)
{
	tmpvars_t *ret = (tmpvars_t *)malloc(sizeof(tmpvars_t));
	ret->blen = blen;
	ret->readlen = readlen;
	ret->bs_ptr = bs_ptr;
	ret->nuc_indices[0] = ret->nuc_indices[1] = 0;
	ret->buffers = (tmpbuffers_t *)malloc(sizeof(tmpbuffers_t));
	ret->buffers->name_buffer[0] = '@';
	ret->buffers->name_buffer[blen] = '\0';
	return ret;
}


static tmpvars_t init_tmpvars(char *bs_ptr, int blen, int readlen)
{
	tmpvars_t ret = {
			.blen = blen,
			.readlen = readlen,
			.bs_ptr = bs_ptr,
			.nuc_indices = {0, 0}
	};
	return ret;
}

static inline void tmpvars_destroy(tmpvars_t *tmp)
{
	free(tmp->buffers);
	free(tmp);
	tmp = NULL;
	return;
}


/*
 * :param: seq - [arg/kseq_t *] a filled-in kseq object.
 * :param: buf - a pre-allocated buffer or malloc'd char_ptr with enough space for the barcode and the null terminus.
 * :returns:
 */
static inline void cp_bs2buf(kseq_t *seq, char *buf)
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



static inline splitterhash_params_t *init_splitterhash_mss(mss_settings_t *settings_ptr, mark_splitter_t *splitter_ptr)
{
#if DBG
	if(!settings_ptr) {
		fprintf(stderr, "Settings struct null. Abort mission!\n");
		exit(EXIT_FAILURE);
	}
	if(!splitter_ptr) {
		fprintf(stderr, "Splitter struct null. Abort mission!\n");
		exit(EXIT_FAILURE);
	}
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
	ret->n = splitter_ptr->n_handles;
	ret->outfnames_r1 = (char **)malloc(ret->n * sizeof(char *));
	ret->outfnames_r2 = (char **)malloc(ret->n * sizeof(char *));
	ret->infnames_r1 = (char **)malloc(ret->n * sizeof(char *));
	ret->infnames_r2 = (char **)malloc(ret->n * sizeof(char *));
	fprintf(stderr, "About to initialize array entries for splitter.\n");
	for(int i = 0; i < ret->n; ++i) {
		if(!splitter_ptr->fnames_r1[i]) {
			fprintf(stderr, "Input r1 filename with index %i null. Abort!\n", i);
			exit(EXIT_FAILURE);
		}
		if(!splitter_ptr->fnames_r2[i]) {
			fprintf(stderr, "Input r2 filename with index %i null. Abort!\n", i);
			exit(EXIT_FAILURE);
		}
		ret->infnames_r1[i] = splitter_ptr->fnames_r1[i];
		ret->infnames_r2[i] = splitter_ptr->fnames_r2[i]; // Does not allocate memory.  This is freed by mark_splitter_t!
		sprintf(tmp_buffer, "%s.%i.R1.dmp.fastq", settings_ptr->output_basename, i);
		ret->outfnames_r1[i] = strdup(tmp_buffer);
		sprintf(tmp_buffer, "%s.%i.R2.dmp.fastq", settings_ptr->output_basename, i);
		ret->outfnames_r2[i] = strdup(tmp_buffer);
	}
	fprintf(stderr, "Finished initializing splitterhash with size %i and output basename %s.\n", ret->n, settings_ptr->output_basename);
	return ret;
}


static inline splitterhash_params_t *init_splitterhash(mssi_settings_t *settings_ptr, mark_splitter_t *splitter_ptr)
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

static inline void splitterhash_destroy(splitterhash_params_t *params)
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


#endif // KHASH_DMP_CORE_H
