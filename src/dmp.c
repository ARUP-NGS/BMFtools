#include "dmp.h"
//Function definitions


static void print_usage(char *argv[]) {
	fprintf(stderr, "Usage: %s -o <output_path> <input_path> (<optional other input_paths>)\n"
					"Set output_path to \"-\" for stdout, "
					"input_path to \"-\" for stdin.\n", argv[0]);
}

static void print_opt_err(char *argv[], char optarg[]) {
	fprintf(stderr, "Invalid argument %s. See usage.\n", optarg);
	print_usage(argv);
	exit(1);
}

// TODO: Use open_memstream to parallelize DMP and merge without a cat!

int main(int argc, char* argv[]) {
#ifdef BMF_THREADS
	int threads = 0;
#endif
	char *outfname = NULL;
	int c;
	char **infnames = (char **)malloc((argc - 3) * sizeof(char *));
#ifdef BMF_THREADS
	while ((c = getopt(argc, argv, "l:o:t:h")) > -1) {
#endif
	while ((c = getopt(argc, argv, "l:o:t:h")) > -1) {
		switch(c) {
#ifdef BMF_THREADS
			case 't': threads = atoi(optarg); break;
#else
			case 't': fprintf(stderr, "bmftools_dmp was compiled without BMF_THREADS. Invalid option!\n"); return 1;
#endif
			case 'o': outfname = strdup(optarg); break;
			case 'h': print_usage(argv); exit(1);
			default: print_opt_err(argv, optarg);
		}
	}
	if(argc < 4) {
		fprintf(stderr, "Too few arguments. See usage.\n");
		print_usage(argv);
		return 1;
	}
	int fnames_count = 0;
	fprintf(stderr, "Entry at optind %i: %s\n", optind, argv[optind]);
	for(int i = optind; i < argc; i++) {
		infnames[fnames_count] = strdup(argv[optind]);
		fnames_count++;
	}
	int abort = 0;
	int retcode;
	int index;
	char *mode = (fnames_count == 1) ? "w": "a"; // Append to the file in the case of looping
	for(index = 0; index < fnames_count; index++) {
		fprintf(stderr, "About to call bmftools_dmp_wrapper for input %s and output %s.\n", infnames[index], outfname);
		retcode = bmftools_dmp_wrapper(infnames[index], outfname);
		if(retcode) {
			abort = 1;
			break;
		}
	}
#if !NDEBUG
	// DEBUG code goes here.
#endif
	if(abort) {
		fprintf(stderr, "bmftools_dmp_wrapper and bmftools_dmp_core returned a non-zero exit status. (EXIT_FAILURE). Abort!\n");
	}
	for(int i = 0; i < fnames_count; i++) {
		free(infnames[i]);
	}
	free(infnames);
	free(outfname);
	return retcode;
}


static int bmftools_dmp_wrapper(char *input_path, char *output_path) {
	/*
	 * Set output_path to NULL to write to stdout.
	 * Set input_path to "-" or NULL to read from stdin.
	 */
	fprintf(stderr, "Starting bmftools_dmp_wrapper.\n");
	FILE *in_handle;
	FILE *out_handle;
	if(input_path[0] == '-' || !input_path) {
		fprintf(stderr, "Input path set to stdin.\n");
		in_handle = stdin;
	}
	else {
		in_handle = fopen(input_path, "r");
	}
	if(!output_path) out_handle = stdout;
	else {
		out_handle = fopen(output_path, "w");
	}
	gzFile fp = gzdopen(fileno(in_handle), "r");
	kseq_t *seq = kseq_init(fp);
	fprintf(stderr, "Opened file handles, initiated kseq parser.\n");
	int retcode = bmftools_dmp_core(seq, out_handle);
	kseq_destroy(seq);
	return retcode;
}


/*
 * TODO:
 * Expand the writing out of KingFisher_t to include the RC count.
 * Only write this out of the count >= 0.
 */

static int bmftools_dmp_core(kseq_t *seq, FILE *out_handle) {
	tmpbuffers_t tmp;
	int l, readlen;
	int nuc_indices[2];
	l = kseq_read(seq);
	fprintf(stderr, "Starting bmftools_dmp_core.\n");
	if(l < 0) {
		fprintf(stderr, "Could not open input fastq. Abort!\n");
		exit(1);
	}
	readlen = strlen(seq->seq.s);
	fprintf(stderr, "read length (inferred): %i\n", readlen);
	KingFisher_t Holloway = init_kf(readlen);
	/*
	KingFisher_t Holloway = init_kf(readlen);
	*/
	KingFisher_t *Hook = &Holloway;
	char *bs_ptr = barcode_mem_view(seq); //Note: NOT null-terminated at the end of the barcode.
	int blen = infer_barcode_length(bs_ptr);
#if !NDDEBUG
	fprintf(stderr, "Barcode length (inferred): %i\n", blen);
#endif
	char *current_barcode = (char *)calloc(blen + 1, 1);
	memcpy(current_barcode, bs_ptr, blen);
#if !NDEBUG
	fprintf(stderr, "Current barcode: %s.\n", current_barcode);
#endif
	pushback_kseq(Hook, seq, nuc_indices, blen); // Initialize Hook with
	while ((l = kseq_read(seq)) >= 0) {
		bs_ptr = barcode_mem_view(seq);
#if !NDEBUG
		if(!bs_ptr) { // bs_ptr is NULL
			destroy_kf(Hook);
			fprintf(stderr, "Malformed fastq comment field - missing the second delimiter. Abort!\n");
			return 1;
		}
#endif
		if(memcmp(bs_ptr, current_barcode, blen) == 0) {
#if !NDEBUG
			fprintf(stderr, "Same barcode. Continue pushing back.\n");
#endif
			pushback_kseq(Hook, seq, nuc_indices, blen);
		}
		else {
#if !NDEBUG
			fprintf(stderr, "Different barcode. Write out result.\n");
#endif
			dmp_process_write(Hook, out_handle, blen, &tmp);
			clear_kf(Hook); // Reset Holloway
			memcpy(current_barcode, bs_ptr, blen * sizeof(char)); // Update working barcode.
			pushback_kseq(Hook, seq, nuc_indices, blen);
		}
	}
	if(Hook->length) { // If length is not 0
		dmp_process_write(Hook, out_handle, blen, &tmp);
	}
	destroy_kf(Hook);
	return 0;
}
