#include "khash_dmp_core.h"

void print_khashdmp_usage(char *argv[]) {
	fprintf(stderr, "Usage: %s -o <output_filename> <input_filename>.\n", argv[0]);
}

void print_khashdmp_opt_err(char *argv[], char *optarg) {
	fprintf(stderr, "Invalid argument %s. See usage.\n", optarg);
	print_khashdmp_usage(argv);
	exit(1);
}

splitterhash_params_t *init_splitterhash_mss(mss_settings_t *settings_ptr, mark_splitter_t *splitter_ptr)
{
	if(!settings_ptr) {
		fprintf(stderr, "[E:%s] Settings pointer null. Abort!\n", __func__);
		exit(EXIT_FAILURE);
	}
	if(!settings_ptr->output_basename) {
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
	for(int i = 0; i < ret->n; ++i) {
		if(!splitter_ptr->fnames_r1[i]) {
			fprintf(stderr, "[E:%s] Input r1 filename with index %i null. Abort!\n", __func__, i);
			exit(EXIT_FAILURE);
		}
		if(!splitter_ptr->fnames_r2[i]) {
			fprintf(stderr, "[E:%s] Input r2 filename with index %i null. Abort!\n", __func__, i);
			exit(EXIT_FAILURE);
		}
		ret->infnames_r1[i] = splitter_ptr->fnames_r1[i];
		ret->infnames_r2[i] = splitter_ptr->fnames_r2[i]; // Does not allocate memory.  This is freed by mark_splitter_t!
		sprintf(tmp_buffer, "%s.%i.R1.dmp.fastq", settings_ptr->output_basename, i);
		ret->outfnames_r1[i] = strdup(tmp_buffer);
		sprintf(tmp_buffer, "%s.%i.R2.dmp.fastq", settings_ptr->output_basename, i);
		ret->outfnames_r2[i] = strdup(tmp_buffer);
	}
	return ret;
}


splitterhash_params_t *init_splitterhash(mssi_settings_t *settings_ptr, mark_splitter_t *splitter_ptr)
{
	if(!settings_ptr) {
		fprintf(stderr, "[E:%s] Settings pointer null. Abort!\n", __func__);
		exit(EXIT_FAILURE);
	}
	if(!settings_ptr->output_basename) {
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
		sprintf(tmp_buffer, "%s.%i.R1.dmp.fastq", settings_ptr->output_basename, i);
		ret->outfnames_r1[i] = strdup(tmp_buffer);
		sprintf(tmp_buffer, "%s.%i.R2.dmp.fastq", settings_ptr->output_basename, i);
		ret->outfnames_r2[i] = strdup(tmp_buffer);
	}
	return ret;
}

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

tmpvars_t *init_tmpvars_p(char *bs_ptr, int blen, int readlen)
{
	tmpvars_t *ret = (tmpvars_t *)malloc(sizeof(tmpvars_t));
	ret->blen = blen;
	ret->readlen = readlen;
	ret->bs_ptr = bs_ptr;
	ret->buffers = (tmpbuffers_t *)malloc(sizeof(tmpbuffers_t));
	ret->buffers->name_buffer[0] = '@';
	ret->buffers->name_buffer[blen] = '\0';
	return ret;
}



int khash_dmp_main(int argc, char *argv[])
{
	char *outfname = NULL;
	char *infname = NULL;
	int c;
	while ((c = getopt(argc, argv, "ho:")) > -1) {
		switch(c) {
			case 'o': outfname = strdup(optarg);break;
			case 'h': print_khashdmp_usage(argv); return 0;
			default: print_khashdmp_opt_err(argv, optarg);
		}
	}
	if(argc < 4) {
		fprintf(stderr, "[E:%s] Required arguments missing. See usage.\n", __func__);
		print_khashdmp_usage(argv);
		exit(1);
	}
	infname = strdup(argv[optind]);

	khash_dmp_core(infname, outfname);
	free(outfname);
	free(infname);
	return 0;
}


void khash_dmp_core(char *infname, char *outfname)
{
	FILE *in_handle;
	FILE *out_handle;
	if(!outfname && !infname) {
		fprintf(stderr, "[E:%s] in_handle and out_handle are both null. Abort mission!\n", __func__);
		exit(EXIT_FAILURE);
	}
	if(infname[0] == '-' || !infname) in_handle = stdin;
	else {
		in_handle = fopen(infname, "r");
	}
	if(!outfname) out_handle = stdout;
	else {
		out_handle = fopen(outfname, "w");
	}
	if(!in_handle) {
		fprintf(stderr, "[E:%s] Could not open %s for reading. Abort mission!\n", __func__, infname);
		exit(EXIT_FAILURE);
	}
	gzFile fp = gzdopen(fileno(in_handle), "r");
	kseq_t *seq = kseq_init(fp);
	// Initialized kseq
	int l = kseq_read(seq);
	if(l < 0) {
		fprintf(stderr, "[%s]: Could not open fastq file (%s). Abort mission!\n",
				__func__, strcmp(infname, "-") == 0 ? "stdin": infname);
		exit(1);
	}
	char *bs_ptr = barcode_mem_view(seq);
	int blen = infer_barcode_length(bs_ptr);
#if !NDEBUG
	fprintf(stderr, "[D:%s] Barcode length (inferred): %i.\n", __func__, blen);
#endif
	tmpvars_t *tmp = init_tmpvars_p(bs_ptr, blen, seq->seq.l);
	memcpy(tmp->key, bs_ptr, blen);
	tmp->key[blen] = '\0';
	// Start hash table
	hk_t *hash = NULL;
	hk_t *current_entry = (hk_t *)malloc(sizeof(hk_t));
	hk_t *tmp_hk = current_entry; // Save the pointer location for later comparison.
	cp_view2buf(bs_ptr, current_entry->id);
	current_entry->value = init_kfp(tmp->readlen);
	HASH_ADD_STR(hash, id, current_entry);
	pushback_kseq(current_entry->value, seq, tmp->blen);

	uint64_t count = 0;
	while(LIKELY((l = kseq_read(seq)) >= 0)) {
		if(++count % 1000000 == 0)
			fprintf(stderr, "[%s] Number of records processed: %lu.\n", __func__, count);
		cp_view2buf(seq->comment.s + 14, tmp->key);
		HASH_FIND_STR(hash, tmp->key, tmp_hk);
		if(!tmp_hk) {
			tmp_hk = (hk_t *)malloc(sizeof(hk_t));
			tmp_hk->value = init_kfp(tmp->readlen);
			cp_view2buf(seq->comment.s + 14, tmp_hk->id);
			pushback_kseq(tmp_hk->value, seq, tmp->blen);
			HASH_ADD_STR(hash, id, tmp_hk);
		}
		else
			pushback_kseq(tmp_hk->value, seq, tmp->blen);
	}
	fprintf(stderr, "[%s] Loaded all fastq records into memory for meta-analysis. Now writing out to file!\n", __func__);
	HASH_ITER(hh, hash, current_entry, tmp_hk) {
		dmp_process_write(current_entry->value, out_handle, tmp->buffers);
		destroy_kf(current_entry->value);
		HASH_DEL(hash, current_entry);
		free(current_entry);
	}
	fprintf(stderr, "[%s] Cleaning up.\n", __func__);
	gzclose(fp);
	fclose(out_handle);
	kseq_destroy(seq);
	tmpvars_destroy(tmp);
}
