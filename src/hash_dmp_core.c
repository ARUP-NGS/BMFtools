#include "hash_dmp_core.h"

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
    ret->buffers->cons_seq_buffer[readlen] = '\0';
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
	if(argc < 2) {
		fprintf(stderr, "[E:%s] Required arguments missing. See usage.\n", __func__);
		print_khashdmp_usage(argv);
		exit(1);
	}
	infname = strdup(argv[optind]);

	hash_dmp_core(infname, outfname);
	free(outfname);
	free(infname);
	return 0;
}


void hash_dmp_core(char *infname, char *outfname)
{
	if(!outfname && !infname)
		fprintf(stderr, "[E:%s] in_handle and out_handle are both null. Abort mission!\n", __func__),
		exit(EXIT_FAILURE);
	FILE *in_handle = (infname[0] == '-' || !infname) ? stdin: fopen(infname, "r");
	FILE *out_handle = (!outfname || *outfname == '-') ? stdout: fopen(outfname, "w");
	if(!in_handle)
		fprintf(stderr, "[E:%s] Could not open %s for reading. Abort mission!\n", __func__, infname),
		exit(EXIT_FAILURE);
	gzFile fp = gzdopen(fileno(in_handle), "r");
	kseq_t *seq = kseq_init(fp);
	// Initialized kseq
	int l = kseq_read(seq);
	if(l < 0)
		fprintf(stderr, "[%s]: Could not open fastq file (%s). Abort mission!\n",
				__func__, strcmp(infname, "-") == 0 ? "stdin": infname),
		exit(1);
	char *bs_ptr = barcode_mem_view(seq);
	const int blen = infer_barcode_length(bs_ptr);
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
	pushback_kseq(current_entry->value, seq, blen);

	uint64_t count = 0;
	while(LIKELY((l = kseq_read(seq)) >= 0)) {
		if(++count % 1000000 == 0)
			fprintf(stderr, "[%s] Number of records read from '%s': %lu.\n", __func__,
					strcmp("-", infname) == 0 ? "stdin": infname,count);
		cp_view2buf(seq->comment.s + 14, tmp->key);
		HASH_FIND_STR(hash, tmp->key, tmp_hk);
		if(!tmp_hk) {
			tmp_hk = (hk_t *)malloc(sizeof(hk_t));
			tmp_hk->value = init_kfp(tmp->readlen);
			cp_view2buf(seq->comment.s + 14, tmp_hk->id);
			pushback_kseq(tmp_hk->value, seq, blen);
			HASH_ADD_STR(hash, id, tmp_hk);
		}
		else
			pushback_kseq(tmp_hk->value, seq, blen);
	}
	fprintf(stderr, "[%s] Loaded all fastq records into memory for meta-analysis. Now writing out to file ('%s')!\n", __func__, outfname);
	HASH_ITER(hh, hash, current_entry, tmp_hk) {
		dmp_process_write(current_entry->value, out_handle, tmp->buffers, 0);
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

void stranded_hash_dmp_core(char *infname, char *outfname)
{
	FILE *in_handle = (infname[0] == '-' || !infname) ? stdin: fopen(infname, "r");
	FILE *out_handle = (!outfname || *outfname == '-') ? stdout: fopen(outfname, "w");
	if(!in_handle)
		fprintf(stderr, "[E:%s] Could not open %s for reading. Abort mission!\n", __func__, infname),
		exit(EXIT_FAILURE);
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
	hk_t *hfor = NULL, *hrev = NULL; // Hash forward, hash reverse
	hk_t *crev = (hk_t *)malloc(sizeof(hk_t)); // Current reverse, current forward.
	hk_t *cfor = (hk_t *)malloc(sizeof(hk_t));
	hk_t *tmp_hkr = crev, *tmp_hkf = cfor;
	if(*bs_ptr == 'F') {
		cp_view2buf(bs_ptr + 1, cfor->id);
		cfor->value = init_kfp(tmp->readlen);
		HASH_ADD_STR(hfor, id, cfor);
		pushback_kseq(cfor->value, seq, blen);
	}
	else {
		cp_view2buf(bs_ptr + 1, crev->id);
		crev->value = init_kfp(tmp->readlen);
		HASH_ADD_STR(hrev, id, crev);
		pushback_kseq(hrev->value, seq, blen);
	}

	uint64_t count = 0;
#if !NDEBUG
	uint64_t fcount = 0, rcount = 0;
#endif
	while(LIKELY((l = kseq_read(seq)) >= 0)) {
		if(UNLIKELY(++count % 1000000 == 0))
			fprintf(stderr, "[%s:%s] Number of records processed: %lu.\n", __func__,
					*infname == '-' ? "stdin" : infname, count);
		if(*(seq->comment.s + 14) == 'F') {
#if !NDEBUG
			++fcount;
#endif
			cp_view2buf(seq->comment.s + 15, tmp->key); // 15 to skip the first tags AND the first base ('F'/'R')
			HASH_FIND_STR(hfor, tmp->key, tmp_hkf);
			if(!tmp_hkf) {
				tmp_hkf = (hk_t *)malloc(sizeof(hk_t));
				tmp_hkf->value = init_kfp(tmp->readlen);
				cp_view2buf(seq->comment.s + 15, tmp_hkf->id);
				pushback_kseq(tmp_hkf->value, seq, blen);
				HASH_ADD_STR(hfor, id, tmp_hkf);
			}
			else
				pushback_kseq(tmp_hkf->value, seq, blen);
		}
		else {
#if !NDEBUG
			++rcount;
#endif
			cp_view2buf(seq->comment.s + 15, tmp->key);
			HASH_FIND_STR(hrev, tmp->key, tmp_hkr);
			if(!tmp_hkr) {
				tmp_hkr = (hk_t *)malloc(sizeof(hk_t));
				tmp_hkr->value = init_kfp(tmp->readlen);
				cp_view2buf(seq->comment.s + 15, tmp_hkr->id);
				pushback_kseq(tmp_hkr->value, seq, blen);
				HASH_ADD_STR(hrev, id, tmp_hkr);
			}
			else
				pushback_kseq(tmp_hkr->value, seq, blen);
		}
	}
#if !NDEBUG
	fprintf(stderr, "[%s] Number of reverse reads: %lu. Number of forward reads: %lu.\n", __func__, rcount, fcount);
#endif
	fprintf(stderr, "[%s] Loaded all fastq records into memory for meta-analysis. Now writing out to file ('%s')!\n", __func__, outfname);
	// Write out all unmatched in forward and handle all barcodes handled from both strands.
	HASH_ITER(hh, hfor, cfor, tmp_hkf) {
		HASH_FIND_STR(hrev, cfor->id, crev);
		if(!crev) {
			dmp_process_write(cfor->value, out_handle, tmp->buffers, 0); // No reverse strand found. \='{
			destroy_kf(cfor->value);
			HASH_DEL(hfor, cfor);
			free(cfor);
			continue;
		}
		stranded_process_write(cfor->value, crev->value, out_handle, tmp->buffers);
		destroy_kf(cfor->value), destroy_kf(crev->value);
		HASH_DEL(hrev, crev); HASH_DEL(hfor, cfor);
		cond_free(crev); cond_free(cfor);
	}
	// Handle the barcodes in reverse not in forward
	HASH_ITER(hh, hrev, crev, tmp_hkr) {
		dmp_process_write(crev->value, out_handle, tmp->buffers, 1); // Only reverse strand found
		destroy_kf(crev->value);
		HASH_DEL(hrev, crev);
		cond_free(crev);
	}
	fprintf(stderr, "[%s] Cleaning up.\n", __func__);
	gzclose(fp);
	fclose(in_handle), fclose(out_handle);
	kseq_destroy(seq);
	tmpvars_destroy(tmp);
}
