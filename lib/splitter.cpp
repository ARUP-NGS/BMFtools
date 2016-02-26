#include "splitter.h"

void splitterhash_destroy(splitterhash_params_t *params)
{
	/*
	 * Is this freed by splitterhash_params?
	for(int i = 0; i < params->n; ++i) {
		LOG_DEBUG("i: %i.\n", i);
		cond_free(params->outfnames_r1[i]);
		cond_free(params->outfnames_r2[i]);
		cond_free(params->infnames_r1[i]);
		cond_free(params->infnames_r2[i]);
	}
	*/
	cond_free(params->outfnames_r1);
	cond_free(params->outfnames_r2);
	cond_free(params->infnames_r1);
	cond_free(params->infnames_r2);
	cond_free(params);
}

void free_marksplit_settings(marksplit_settings_t settings)
{
	cond_free(settings.tmp_basename);
	cond_free(settings.input_r1_path);
	cond_free(settings.input_r2_path);
	cond_free(settings.index_fq_path);
	cond_free(settings.rescaler);
	cond_free(settings.rescaler_path);
    cond_free(settings.homing_sequence);
    cond_free(settings.ffq_prefix);
}

void free_marksplit_settings_ptr(marksplit_settings_t *settings)
{
	free_marksplit_settings(*settings);
	free(settings);
}

splitterhash_params_t *init_splitterhash(marksplit_settings_t *settings_ptr, mark_splitter_t *splitter_ptr)
{
	if(!settings_ptr) {
		LOG_EXIT("Settings pointer null. Abort!\n");
	}
	if(!settings_ptr->tmp_basename) {
		fprintf(stderr, "[E:%s] Output basename not set. Abort!\n", __func__);
		exit(EXIT_FAILURE);
	}
	if(!splitter_ptr) {
		fprintf(stderr, "[E:%s] Splitter pointer null. Abort!\n", __func__);
		exit(EXIT_FAILURE);
	}
	kstring_t ks = {0, 0, NULL};
	splitterhash_params_t *ret = (splitterhash_params_t *)malloc(sizeof(splitterhash_params_t));
	ret->n = splitter_ptr->n_handles;
	if(settings_ptr->is_se) {
		ret->outfnames_r1 = (char **)malloc(ret->n * sizeof(char *));
		ret->infnames_r1 = (char **)malloc(ret->n * sizeof(char *));
		for(int i = 0; i < splitter_ptr->n_handles; ++i) {
			ret->infnames_r1[i] = splitter_ptr->fnames_r1[i];
			ks.l = 0;
			ksprintf(&ks, "%s.%i.dmp.fastq", settings_ptr->tmp_basename, i);
			ret->outfnames_r1[i] = ks.s;
		}
	} else {
		ret->outfnames_r1 = (char **)malloc(ret->n * sizeof(char *));
		ret->outfnames_r2 = (char **)malloc(ret->n * sizeof(char *));
		ret->infnames_r1 = (char **)malloc(ret->n * sizeof(char *));
		ret->infnames_r2 = (char **)malloc(ret->n * sizeof(char *));
		for(int i = 0; i < splitter_ptr->n_handles; ++i) {
			ret->infnames_r1[i] = splitter_ptr->fnames_r1[i];
			ret->infnames_r2[i] = splitter_ptr->fnames_r2[i]; // Does not allocate memory.  This is freed by mark_splitter_t!
			ks.l = 0;
			ksprintf(&ks, "%s.%i.R1.dmp.fastq", settings_ptr->tmp_basename, i);
			ret->outfnames_r1[i] = strdup(ks.s);
			ks.l = 0;
			ksprintf(&ks, "%s.%i.R2.dmp.fastq", settings_ptr->tmp_basename, i);
			ret->outfnames_r2[i] = strdup(ks.s);
		}
	}
	free(ks.s);
	return ret;
}


void splitter_destroy(mark_splitter_t *var)
{
	for(int i = 0; i < var->n_handles; i++) {
		cond_free(var->fnames_r1[i]); cond_free(var->fnames_r2[i]);
	}
	cond_free(var->tmp_out_handles_r1);
	cond_free(var->tmp_out_handles_r2);
	cond_free(var);
}


mark_splitter_t init_splitter_pe(marksplit_settings_t* settings_ptr)
{
	mark_splitter_t ret = {
		NULL, // tmp_out_handles_r1
		NULL, // tmp_out_handles_r2
		settings_ptr->n_nucs, // n_nucs
		(int)ipow(4, settings_ptr->n_nucs), // n_handles
		NULL, // infnames_r1
		NULL  // infnames_r2
	};
	ret.tmp_out_handles_r1 = (gzFile *)malloc(settings_ptr->n_handles * sizeof(gzFile));
	ret.tmp_out_handles_r2 = (gzFile *)malloc(settings_ptr->n_handles * sizeof(gzFile));
	ret.fnames_r1 = (char **)malloc(ret.n_handles * sizeof(char *));
	ret.fnames_r2 = (char **)malloc(ret.n_handles * sizeof(char *));
	kstring_t ks = {0, 0, NULL};
	for (int i = 0; i < ret.n_handles; i++) {
		ks.l = 0;
		ksprintf(&ks, "%s.tmp.%i.R1.fastq", settings_ptr->tmp_basename, i);
		ret.fnames_r1[i] = kstrdup(&ks);
		ks.l = 0;
		ksprintf(&ks, "%s.tmp.%i.R2.fastq", settings_ptr->tmp_basename, i);
		ret.fnames_r2[i] = kstrdup(&ks);
		ret.tmp_out_handles_r1[i] = gzopen(ret.fnames_r1[i], settings_ptr->mode);
		ret.tmp_out_handles_r2[i] = gzopen(ret.fnames_r2[i], settings_ptr->mode);
	}
	return ret;
}

mark_splitter_t init_splitter_se(marksplit_settings_t* settings_ptr)
{
	mark_splitter_t ret = {
        NULL, // tmp_out_handles_r1
        NULL, // tmp_out_handles_r2
        settings_ptr->n_nucs, // n_nucs
		(int)ipow(4, settings_ptr->n_nucs), // n_handles
        NULL, // infnames_r1
        NULL  // infnames_r2
	};
	ret.tmp_out_handles_r1 = (gzFile *)malloc(settings_ptr->n_handles * sizeof(gzFile));
	ret.fnames_r1 = (char **)malloc(ret.n_handles * sizeof(char *));
	kstring_t ks = {0, 0, NULL};
	for (int i = 0; i < ret.n_handles; i++) {
		ks.l = 0;
		ksprintf(&ks, "%s.tmp.%i.fastq", settings_ptr->tmp_basename, i);
		ret.fnames_r1[i] = kstrdup(&ks);
		ret.tmp_out_handles_r1[i] = gzopen(ret.fnames_r1[i], settings_ptr->mode);
	}
	return ret;
}

mark_splitter_t init_splitter(marksplit_settings_t* settings_ptr)
{
	if(settings_ptr->is_se) {
#if !NDEBUG
		fprintf(stderr, "[D:%s] Initializing single-end splitter.\n", __func__);
#endif
		return init_splitter_se(settings_ptr);
	} else {
#if !NDEBUG
		fprintf(stderr, "[D:%s] Initializing paired-end splitter.\n", __func__);
#endif
		return init_splitter_pe(settings_ptr);
	}
}

#undef INIT_SPLITTER