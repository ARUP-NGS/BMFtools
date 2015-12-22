#include "splitter.h"

void splitterhash_destroy(splitterhash_params_t *params)
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

void free_marksplit_settings(marksplit_settings_t settings)
{
	cond_free(settings.tmp_basename);
	cond_free(settings.input_r1_path);
	cond_free(settings.input_r2_path);
	cond_free(settings.index_fq_path);
	cond_free(settings.rescaler);
	cond_free(settings.rescaler_path);
}

void free_marksplit_settings_ptr(marksplit_settings_t *settings)
{
	free_marksplit_settings(*settings);
	free(settings);
}

splitterhash_params_t *init_splitterhash(marksplit_settings_t *settings_ptr, mark_splitter_t *splitter_ptr)
{
	if(!settings_ptr) {
		fprintf(stderr, "[E:%s] Settings pointer null. Abort!\n", __func__);
		exit(EXIT_FAILURE);
	}
	if(!settings_ptr->tmp_basename) {
		fprintf(stderr, "[E:%s] Output basename not set. Abort!\n", __func__);
		exit(EXIT_FAILURE);
	}
	if(!splitter_ptr) {
		fprintf(stderr, "[E:%s] Splitter pointer null. Abort!\n", __func__);
		exit(EXIT_FAILURE);
	}
	char tmp_buffer [METASYNTACTIC_FNAME_BUFLEN];
	splitterhash_params_t *ret = (splitterhash_params_t *)malloc(sizeof(splitterhash_params_t));
	ret->n = splitter_ptr->n_handles;
	ret->outfnames_r1 = (char **)malloc(ret->n * sizeof(char *));
	ret->outfnames_r2 = (char **)malloc(ret->n * sizeof(char *));
	ret->infnames_r1 = (char **)malloc(ret->n * sizeof(char *));
	ret->infnames_r2 = (char **)malloc(ret->n * sizeof(char *));
	for(int i = 0; i < splitter_ptr->n_handles; ++i) {
		ret->infnames_r1[i] = splitter_ptr->fnames_r1[i];
		ret->infnames_r2[i] = splitter_ptr->fnames_r2[i]; // Does not allocate memory.  This is freed by mark_splitter_t!
		sprintf(tmp_buffer, "%s.%i.R1.dmp.fastq", settings_ptr->tmp_basename, i);
		ret->outfnames_r1[i] = strdup(tmp_buffer);
		sprintf(tmp_buffer, "%s.%i.R2.dmp.fastq", settings_ptr->tmp_basename, i);
		ret->outfnames_r2[i] = strdup(tmp_buffer);
	}
	return ret;
}


void splitter_destroy(mark_splitter_t *var)
{
	for(int i = 0; i < var->n_handles; i++) {
		cond_free(var->fnames_r1[i]); cond_free(var->fnames_r2[i]);
	}
	free(var->tmp_out_handles_r1);
	free(var->tmp_out_handles_r2);
	free(var);
}

#define INIT_SPLITTER(settings_ptr) \
	do {\
		mark_splitter_t ret = {\
			.n_handles = ipow(4, settings_ptr->n_nucs),\
			.n_nucs = settings_ptr->n_nucs,\
			.fnames_r1 = NULL,\
			.fnames_r2 = NULL,\
			.tmp_out_handles_r1 = NULL,\
			.tmp_out_handles_r2 = NULL\
		};\
		ret.tmp_out_handles_r1 = (FILE **)malloc(settings_ptr->n_handles * sizeof(FILE *));\
		ret.tmp_out_handles_r2 = (FILE **)malloc(settings_ptr->n_handles * sizeof(FILE *));\
		ret.fnames_r1 = (char **)malloc(ret.n_handles * sizeof(char *));\
		ret.fnames_r2 = (char **)malloc(ret.n_handles * sizeof(char *));\
		char tmp_buffer [METASYNTACTIC_FNAME_BUFLEN];\
		for (int i = 0; i < ret.n_handles; i++) {\
			sprintf(tmp_buffer, "%s.tmp.%i.R1.fastq", settings_ptr->tmp_basename, i);\
			ret.fnames_r1[i] = strdup(tmp_buffer);\
			sprintf(tmp_buffer, "%s.tmp.%i.R2.fastq", settings_ptr->tmp_basename, i);\
			ret.fnames_r2[i] = strdup(tmp_buffer);\
			ret.tmp_out_handles_r1[i] = fopen(ret.fnames_r1[i], "w");\
			ret.tmp_out_handles_r2[i] = fopen(ret.fnames_r2[i], "w");\
		}\
		return ret;\
	} while(0)

mark_splitter_t init_splitter(marksplit_settings_t* settings_ptr)
{
	INIT_SPLITTER(settings_ptr);
}

#undef INIT_SPLITTER
