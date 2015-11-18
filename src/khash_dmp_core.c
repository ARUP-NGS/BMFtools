#include "khash_dmp_core.h"

static inline void print_khashdmp_usage(char *argv[]) {
	fprintf(stderr, "Usage: %s -o <output_filename> <input_filename>.\n", argv[0]);
}

static inline void print_khashdmp_opt_err(char *argv[], char *optarg) {
	fprintf(stderr, "Invalid argument %s. See usage.\n", optarg);
	print_khashdmp_usage(argv);
	exit(1);
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
#if DBG
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
#if DBG
	fprintf(stderr, "[khash_dmp_core]: Barcode length (inferred): %i.\n", blen);
#endif
	tmpvars_t *tmp = init_tmpvars_p(bs_ptr, blen, seq->seq.l);
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
	while((l = kseq_read(seq)) >= 0) {
		/*
#if DBG
		fprintf(stderr,"Barcode sequence: %s. Comment: %s.", tmp->bs_ptr, seq->comment.s);
		exit(EXIT_SUCCESS);
		if(!tmp->bs_ptr) {
			fprintf(stderr, "Couldn't get barcode mem view to work for read %s.\n", seq->name.s);
			exit(EXIT_FAILURE);
		}
		fprintf(stderr, "Seq name: %s. Comment: %s. View: %s.\n", seq->name.s, seq->comment.s, tmp->bs_ptr);
		if(!seq->seq.s) {
			fprintf(stderr, "The sequence is somehow null. Stopping reading the file now. \n");
			break;
		}
#endif
		*/
		cp_view2buf(seq->comment.s + 14, tmp->key);
		ki = kh_get(dmp, hash, tmp->key);
		if(ki == kh_end(hash)) {
			ki = kh_put(dmp, hash, tmp->key, &khr);
			kh_val(hash, ki) = init_kfp(tmp->readlen);
		}
#if !NDEBUG
		fprintf(stderr, "Seq: %s. Comment: %s. Quality: %s.\n", seq->seq.s, seq->comment.s, seq->qual.s);
#endif
		pushback_kseq(kh_val(hash, ki), seq, tmp->nuc_indices, tmp->blen);
	}
	fprintf(stderr, "[khash_dmp_core]: Loaded all fastq records into memory for meta-analysis. Now writing out to file!\n");
	for(ki = kh_begin(hash); ki != kh_end(hash); ++ki) {
		if(!kh_exist(hash, ki))
			continue;
		//strcpy(tmp->key, kh_val(hash, ki)->barcode);
		dmp_process_write(kh_val(hash, ki), out_handle, tmp);
		destroy_kf(kh_val(hash, ki));
		kh_val(hash, ki) = NULL;
	}
	kh_destroy(dmp, hash);
	gzclose(fp);
	fclose(out_handle);
	kseq_destroy(seq);
	tmpvars_destroy(tmp);
}
