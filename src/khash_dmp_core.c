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


splitterhash_params_t *init_splitterhash(mssi_settings_t *settings_ptr, mark_splitter_t *splitter_ptr)
{
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
		fprintf(stderr, "Required arguments missing. See usage.\n");
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
		fprintf(stderr, "in_handle and out_handle are both null. Abort mission!\n");
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
		fprintf(stderr, "Could not open %s for reading. Abort mission!\n", infname);
	}
	gzFile fp = gzdopen(fileno(in_handle), "r");
	kseq_t *seq = kseq_init(fp);
#if !NDEBUG
	fprintf(stderr, "[khash_dmp_core]: Opened file handles, initiated kseq parser.\n");
#endif
	// Initialized kseq
	int l = kseq_read(seq);
	if(l < 0) {
		fprintf(stderr, "[khash_dmp_core]: Could not open fastq file (%s). Abort mission!\n",
				strcmp(infname, "-") == 0 ? "stdin": infname);
		exit(1);
	}
	char *bs_ptr = barcode_mem_view(seq);
	int blen = infer_barcode_length(bs_ptr);
#if !NDEBUG
	fprintf(stderr, "[khash_dmp_core]: Barcode length (inferred): %i.\n", blen);
#endif
	tmpvars_t *tmp = init_tmpvars_p(bs_ptr, blen, seq->seq.l);
	tmpbuffers_t *bufs = tmp->buffers;
	memcpy(tmp->key, bs_ptr, blen);
	tmp->key[blen] = '\0';
	// Start hash table
	khash_t(dmp) *hash = kh_init(dmp);
	kh_resize(dmp, hash, 1<<18);
	khiter_t ki;
	int khr;
	ki = kh_put(dmp, hash, tmp->key, &khr);
	kh_val(hash, ki) = init_kfp(tmp->readlen);
	pushback_kseq(kh_val(hash, ki), seq, tmp->nuc_indices, tmp->blen);
	while(LIKELY((l = kseq_read(seq)) >= 0)) {
		cp_view2buf(seq->comment.s + 14, tmp->key);
		if((ki = kh_get(dmp, hash, tmp->key)) == kh_end(hash)) {
			ki = kh_put(dmp, hash, tmp->key, &khr);
			kh_val(hash, ki) = init_kfp(tmp->readlen);
		}
		pushback_kseq(kh_val(hash, ki), seq, tmp->nuc_indices, tmp->blen);
	}
	fprintf(stderr, "[khash_dmp_core]: Loaded all fastq records into memory for meta-analysis. Now writing out to file!\n");
	for(ki = kh_begin(hash); ki != kh_end(hash); ++ki) {
		if(!kh_exist(hash, ki))
			continue;
		//strcpy(tmp->key, kh_val(hash, ki)->barcode);
		dmp_process_write(kh_val(hash, ki), out_handle, bufs);
		destroy_kf(kh_val(hash, ki));
		kh_val(hash, ki) = NULL;
	}
	kh_destroy(dmp, hash);
	gzclose(fp);
	fclose(out_handle);
	kseq_destroy(seq);
	tmpvars_destroy(tmp);
}
