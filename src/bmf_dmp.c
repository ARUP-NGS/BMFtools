/*
 * Conditional reverse complement Rescale Mark Split
 * For inline Loeb adapters only.
 */

#include "bmf_dmp.h"


void print_crms_usage(char *executable)
{
		fprintf(stderr, "Usage: %s <options> <Fq.R1.seq> <Fq.R2.seq>"
						"\nFlags:\n"
						"-l: Number of nucleotides at the beginning of each read to "
						"use for barcode. Final barcode length is twice this. REQUIRED.\n"
						"-s: homing sequence. REQUIRED.\n"
						"-o: Mark/split output basename. Defaults to a random string.\n"
						"-t: Homopolymer failure threshold. A molecular barcode with"
						" a homopolymer of length >= this limit is flagged as QC fail."
						"Default: 10.\n"
						"-n: Number of nucleotides at the beginning of the barcode to use to split the output. Default: 4.\n"
						"-m: Mask first n nucleotides in read for barcode. Default: 1.\n"
						"-p: Number of threads to use if running uthash_dmp.\n"
						"-d: Use this flag to to run hash_dmp.\n"
						"-f: If running hash_dmp, this sets the Final Fastq Prefix. \n"
						"The Final Fastq files will be named '<ffq_prefix>.R1.fq' and '<ffq_prefix>.R2.fq'.\n"
						"-r: Path to flat text file with rescaled quality scores. If not provided, it will not be used.\n"
						"-v: Maximum barcode length for a variable length barcode dataset. If left as default value,"
						" (-1), other barcode lengths will not be considered.\n"
						"-z: Flag to optionally pipe to gzip while producing final fastqs. Default: False.\n"
						"-g: Gzip compression ratio if piping to gzip (-z). Default: 1 (mostly to reduce I/O).\n"
						"-c: Flag to optionally cat all files together in one command. Faster than sequential cats, but might break."
						"In addition, won't work for enormous filenames or too many arguments. Default: False.\n"
						"-u: Set notification/update interval for split. Default: 1000000.\n"
						"-w: Set flag to leave temporary files. Primarily for debugging.\n"
						"-h: Print usage.\n", executable);
}

void print_crms_opt_err(char *arg, char *optarg, char optopt)
{
	print_crms_usage(arg);
	fprintf(stderr, "[E:%s] Unrecognized option %s for flag %c. Abort!\n", __func__, optarg, optopt);
	exit(EXIT_FAILURE);
}

/*
 * Pre-processes (pp) and splits fastqs with inline barcodes.
 */
mark_splitter_t *pp_split_annealed(marksplit_settings_t *settings)
{
#if WRITE_BARCODE_FQ
	FILE *fp = fopen("tmp.molbc.fq", "w");
#endif
	fprintf(stderr, "[%s] Opening fastqs %s and %s.\n", __func__, settings->input_r1_path,
            settings->input_r2_path);
	if(!(strcmp(settings->input_r1_path, settings->input_r2_path))) {
		fprintf(stderr, "Hey, it looks like you're trying to use the same path for both r1 and r2. At least try to fool me by making a symbolic link.\n");
		exit(EXIT_FAILURE);
	}
	if(settings->rescaler_path) {
		settings->rescaler = parse_1d_rescaler(settings->rescaler_path);
	}
	mark_splitter_t *splitter = (mark_splitter_t *)malloc(sizeof(mark_splitter_t));
	*splitter = init_splitter(settings);
	gzFile fp1 = gzopen(settings->input_r1_path, "r");
	gzFile fp2 = gzopen(settings->input_r2_path, "r");
	kseq_t *seq1 = kseq_init(fp1);
	kseq_t *seq2 = kseq_init(fp2);
	mseq_t *rseq1, *rseq2;
	int l1, l2;
	int count = 0;
	char pass_fail;
	settings->blen1_2 = settings->blen / 2;
	if((l1 = kseq_read(seq1)) < 0 || (l2 = kseq_read(seq2)) < 0) {
			fprintf(stderr, "[E:%s] Could not open fastqs for reading. Abort!\n", __func__);
            free_marksplit_settings(*settings);
			splitter_destroy(splitter);
			exit(EXIT_FAILURE);
	}
	fprintf(stderr, "[%s]Read length (inferred): %lu.\n", __func__, seq1->seq.l);
#if !NDEBUG
	int arr_size = seq1->seq.l * 4 * 2 * nqscores;
	if(settings->rescaler) {
		for(int i = 0; i < arr_size; ++i) {
			if(settings->rescaler[i] < 0) {
				fprintf(stderr, "[D:%s] Rescaler's got a negative number"
                        " in pp_split_inline. WTF? %i. Index: %i.\n", __func__, settings->rescaler[i], i);
				exit(EXIT_FAILURE);
			}
			else if(settings->rescaler[i] == 0) {
				fprintf(stderr, "[D:%s] Rescaler's got a zero"
                        " in pp_split_inline. WTF? %i. Index: %i.\n", __func__,settings->rescaler[i], i);
				exit(EXIT_FAILURE);
			}
		}
	}
#endif
	tmp_mseq_t *tmp = init_tm_ptr(seq1->seq.l, settings->blen);
	int default_nlen = settings->blen1_2 + settings->offset + settings->homing_sequence_length;
	int n_len = nlen_homing_default(seq1, seq2, settings, default_nlen, &pass_fail);
	rseq1 = mseq_rescale_init(seq1, settings->rescaler, tmp, 0);
	rseq2 = mseq_rescale_init(seq2, settings->rescaler, tmp, 1);
	rseq1->barcode[settings->blen] = '\0';
	if(pass_fail == '1' && (!test_hp_inline(rseq1->barcode, settings->blen, settings->hp_threshold)))
		pass_fail = '0';
	memcpy(rseq1->barcode + settings->blen1_2, seq2->seq.s + settings->offset, settings->blen1_2);
	uint64_t bin = get_binner_type(rseq1->barcode, settings->n_nucs, uint64_t);
	mask_mseq(rseq1, n_len);
	mask_mseq(rseq2, n_len);
	// Get first barcode.
	update_mseq(rseq1, seq1, settings->rescaler, tmp, n_len, 0, 0);
	update_mseq(rseq2, seq2, settings->rescaler, tmp, n_len, 1, 0);
	mseq2fq_inline(splitter->tmp_out_handles_r1[bin], rseq1, pass_fail, rseq1->barcode);
	mseq2fq_inline(splitter->tmp_out_handles_r2[bin], rseq2, pass_fail, rseq1->barcode);
	do {
		if(UNLIKELY(++count % settings->notification_interval == 0))
			fprintf(stderr, "[%s] Number of records processed: %i.\n", __func__, count);
		// Iterate through second fastq file.
		n_len = nlen_homing_default(seq1, seq2, settings, default_nlen, &pass_fail);
		//fprintf(stdout, "Randomly testing to see if the reading is working. %s", seq1->seq.s);
		update_mseq(rseq1, seq1, settings->rescaler, tmp, n_len, 0, 0);
		update_mseq(rseq2, seq2, settings->rescaler, tmp, n_len, 1, 0);
		memcpy(rseq1->barcode, seq1->seq.s + settings->offset, settings->blen1_2);
		memcpy(rseq1->barcode + settings->blen1_2, seq2->seq.s + settings->offset, settings->blen1_2);
	    if(pass_fail == '1' && (!test_hp_inline(rseq1->barcode, settings->blen, settings->hp_threshold)))
			pass_fail = '0';
		bin = get_binner_type(rseq1->barcode, settings->n_nucs, uint64_t);
		mseq2fq_inline(splitter->tmp_out_handles_r1[bin], rseq1, pass_fail, rseq1->barcode);
		mseq2fq_inline(splitter->tmp_out_handles_r2[bin], rseq2, pass_fail, rseq1->barcode);
	} while (LIKELY(LIKELY(((l1 = kseq_read(seq1)) >= 0)) && LIKELY((l2 = kseq_read(seq2)) >= 0)));
	for(int i = 0; i < splitter->n_handles; ++i) {
		fclose(splitter->tmp_out_handles_r1[i]);
		fclose(splitter->tmp_out_handles_r2[i]);
	}
	tm_destroy(tmp);
	mseq_destroy(rseq1), mseq_destroy(rseq2);
	kseq_destroy(seq1), kseq_destroy(seq2);
	gzclose(fp1), gzclose(fp2);
#if WRITE_BARCODE_FQ
	if(fclose(fp))
		fprintf(stderr, "[E:%s] Failed to close file for barcode fastq.\n"),
		exit(EXIT_FAILURE);
#endif
	return splitter;
}


/*
 * Pre-processes (pp) and splits fastqs with inline barcodes.
 */
mark_splitter_t *pp_split_inline(marksplit_settings_t *settings)
{
#if WRITE_BARCODE_FQ
	FILE *fp = fopen("tmp.molbc.fq", "w");
#endif
	fprintf(stderr, "[%s] Opening fastq files %s and %s.\n", __func__, settings->input_r1_path, settings->input_r2_path);
	if(!(strcmp(settings->input_r1_path, settings->input_r2_path))) {
		fprintf(stderr, "[E:%s] Hey, it looks like you're trying to use the same path for both r1 and r2. "
                "At least try to fool me by making a symbolic link.\n",
                __func__);
		exit(EXIT_FAILURE);
	}
	if(!isfile(settings->input_r1_path) || !isfile(settings->input_r2_path))
		fprintf(stderr, "[E:%s] Could not open read paths: at least one is not a file.\n", __func__),
		exit(EXIT_FAILURE);
	if(settings->rescaler_path) {
		settings->rescaler = parse_1d_rescaler(settings->rescaler_path);
	}
	mark_splitter_t *splitter = (mark_splitter_t *)malloc(sizeof(mark_splitter_t));
	*splitter = init_splitter(settings);
	gzFile fp1 = gzopen(settings->input_r1_path, "r");
	gzFile fp2 = gzopen(settings->input_r2_path, "r");
	kseq_t *seq1 = kseq_init(fp1);
	kseq_t *seq2 = kseq_init(fp2);
	mseq_t *rseq1, *rseq2;
	int l1, l2;
	int count = 0;
	char pass_fail;
	settings->blen1_2 = settings->blen / 2;
	if((l1 = kseq_read(seq1)) < 0 || (l2 = kseq_read(seq2)) < 0) {
			fprintf(stderr, "[E:%s] Could not open fastqs for reading. Abort!\n", __func__);
            free_marksplit_settings(*settings);
            splitter_destroy(splitter);
			exit(EXIT_FAILURE);
	}
	fprintf(stderr, "[%s] Read length (inferred): %lu.\n", __func__, seq1->seq.l);
#if !NDEBUG
	int arr_size = seq1->seq.l * 4 * 2 * nqscores;
	if(settings->rescaler) {
		for(int i = 0; i < arr_size; ++i) {
			if(settings->rescaler[i] < 0) {
				fprintf(stderr, "[E:%s] Rescaler's got a negative number in pp_split_inline."
                        " WTF? %i. Index: %i.\n", __func__, settings->rescaler[i], i);
				exit(EXIT_FAILURE);
			}
			else if(settings->rescaler[i] == 0) {
				fprintf(stderr, "[E:%s] Rescaler's got a zero in pp_split_inline. WTF? %i. Index: %i.\n",
                        __func__, settings->rescaler[i], i);
				exit(EXIT_FAILURE);
			}
		}
	}
#endif
	tmp_mseq_t *tmp = init_tm_ptr(seq1->seq.l, settings->blen);
	int switch_reads = switch_test(seq1, seq2, settings->offset);
	int default_nlen = settings->blen1_2 + settings->offset + settings->homing_sequence_length;
	int n_len = nlen_homing_default(seq1, seq2, settings, default_nlen, &pass_fail);
	rseq1 = mseq_rescale_init(seq1, settings->rescaler, tmp, 0);
	rseq2 = mseq_rescale_init(seq2, settings->rescaler, tmp, 1);
	rseq1->barcode[settings->blen] = '\0';
	if(pass_fail == '1' && (!test_hp_inline(rseq1->barcode, settings->blen, settings->hp_threshold)))
		pass_fail = '0';
	if(switch_reads) {
		memcpy(rseq1->barcode, seq2->seq.s + settings->offset, settings->blen1_2);
		memcpy(rseq1->barcode + settings->blen1_2, seq1->seq.s + settings->offset, settings->blen1_2);
	}
	else {
		memcpy(rseq1->barcode, seq1->seq.s + settings->offset, settings->blen1_2);
		memcpy(rseq1->barcode + settings->blen1_2, seq2->seq.s + settings->offset, settings->blen1_2);
	}
	mask_mseq(rseq1, n_len);
	mask_mseq(rseq2, n_len);
	// Get first barcode.
	update_mseq(rseq1, seq1, settings->rescaler, tmp, n_len, 0, switch_reads);
	update_mseq(rseq2, seq2, settings->rescaler, tmp, n_len, 1, switch_reads);
	uint64_t bin = get_binner_type(rseq1->barcode, settings->n_nucs, uint64_t);
	if(switch_reads) {
		mseq2fq_stranded(splitter->tmp_out_handles_r1[bin], rseq2, pass_fail, rseq1->barcode, 'R');
		mseq2fq_stranded(splitter->tmp_out_handles_r2[bin], rseq1, pass_fail, rseq1->barcode, 'R');
	}
	else {
		mseq2fq_stranded(splitter->tmp_out_handles_r1[bin], rseq1, pass_fail, rseq1->barcode, 'F');
		mseq2fq_stranded(splitter->tmp_out_handles_r2[bin], rseq2, pass_fail, rseq1->barcode, 'F');
	}
	do {
		if(UNLIKELY(++count % settings->notification_interval == 0))
			fprintf(stderr, "[%s] Number of records processed: %i.\n", __func__, count);
		// Iterate through second fastq file.
		n_len = nlen_homing_default(seq1, seq2, settings, default_nlen, &pass_fail);
		switch_reads = switch_test(seq1, seq2, settings->offset);
		//fprintf(stdout, "Randomly testing to see if the reading is working. %s", seq1->seq.s);
		update_mseq(rseq1, seq1, settings->rescaler, tmp, n_len, 0, switch_reads);
		update_mseq(rseq2, seq2, settings->rescaler, tmp, n_len, 1, switch_reads);
		if(switch_reads) {
			memcpy(rseq1->barcode, seq2->seq.s + settings->offset, settings->blen1_2);
			memcpy(rseq1->barcode + settings->blen1_2, seq1->seq.s + settings->offset, settings->blen1_2);
		}
		else {
			memcpy(rseq1->barcode, seq1->seq.s + settings->offset, settings->blen1_2);
			memcpy(rseq1->barcode + settings->blen1_2, seq2->seq.s + settings->offset, settings->blen1_2);
		}
		if(pass_fail - '0' && (!test_hp_inline(rseq1->barcode, settings->blen, settings->hp_threshold)))
			pass_fail = '0';
	    bin = get_binner_type(rseq1->barcode, settings->n_nucs, uint64_t);
		if(switch_reads) {
			mseq2fq_stranded(splitter->tmp_out_handles_r1[bin], rseq2, pass_fail, rseq1->barcode, 'R');
			mseq2fq_stranded(splitter->tmp_out_handles_r2[bin], rseq1, pass_fail, rseq1->barcode, 'R');
		}
		else {
			mseq2fq_stranded(splitter->tmp_out_handles_r1[bin], rseq1, pass_fail, rseq1->barcode, 'F');
			mseq2fq_stranded(splitter->tmp_out_handles_r2[bin], rseq2, pass_fail, rseq1->barcode, 'F');
		}
	} while (LIKELY(LIKELY((l1 = kseq_read(seq1)) >= 0) && LIKELY((l2 = kseq_read(seq2)) >= 0)));
    fprintf(stderr, "[%s] Cleaning up.\n", __func__);
	for(int i = 0; i < splitter->n_handles; ++i) {
		fclose(splitter->tmp_out_handles_r1[i]);
		fclose(splitter->tmp_out_handles_r2[i]);
	}
	tm_destroy(tmp);
	mseq_destroy(rseq1), mseq_destroy(rseq2);
	kseq_destroy(seq1), kseq_destroy(seq2);
	gzclose(fp1), gzclose(fp2);
	return splitter;
}


int crms_main(int argc, char *argv[])
{
    if(argc == 1) print_crms_usage(argv[0]), exit(EXIT_FAILURE);
	// Build settings struct
	marksplit_settings_t settings = {
		.hp_threshold = 10,
		.n_nucs = 2,
		.output_basename = NULL,
		.input_r1_path = NULL,
		.input_r2_path = NULL,
		.index_fq_path = NULL, // This is unused for inline experiments.
		.homing_sequence = NULL,
		.n_handles = 0,
		.notification_interval = 1000000,
		.blen = 0,
		.homing_sequence_length = 0,
		.offset = 1,
		.rescaler = NULL,
		.rescaler_path = NULL,
		.run_hash_dmp = 0,
		.ffq_prefix = NULL,
		.threads = 1,
		.max_blen = -1,
		.gzip_output = 0,
		.panthera = 0,
		.gzip_compression = 1,
		.cleanup = 1,
		.annealed = 0,
		.salt = 0 // This is unused for inline experiments
	};
	//omp_set_dynamic(0); // Tell omp that I want to set my number of threads 4realz
	int c;
	while ((c = getopt(argc, argv, "t:o:n:s:l:m:r:p:f:v:u:g:i:zwcdh?")) > -1) {
		switch(c) {
			case 'c': settings.panthera = 1; break;
			case 'd': settings.run_hash_dmp = 1; break;
			case 'f': settings.ffq_prefix = strdup(optarg); break;
			case 'g': settings.gzip_compression = atoi(optarg); break;
			case '?': // Fall-through
			case 'h': print_crms_usage(argv[0]); return EXIT_SUCCESS;
			case 'l': settings.blen = 2 * atoi(optarg); break;
			case 'm': settings.offset = atoi(optarg); break;
			case 'n': settings.n_nucs = atoi(optarg); break;
			case 'o': settings.output_basename = strdup(optarg); break;
			case 'p': settings.threads = atoi(optarg); break;
			case 'r': settings.rescaler_path = strdup(optarg); break;
			case 's': settings.homing_sequence = strdup(optarg); settings.homing_sequence_length = strlen(settings.homing_sequence); break;
			case 't': settings.hp_threshold = atoi(optarg); break;
			case 'u': settings.notification_interval = atoi(optarg); break;
			case 'v': settings.max_blen = atoi(optarg); break;
			case 'w': settings.cleanup = 0;
			case 'z': settings.gzip_output = 1; break;
			default: print_crms_opt_err(argv[0], optarg, optopt);
		}
	}

	increase_nofile_limit(settings.threads);
	omp_set_num_threads(settings.threads);

	if(argc < 5)
		print_crms_usage(argv[0]), exit(EXIT_FAILURE);

	settings.n_handles = ipow(4, settings.n_nucs);
	if(settings.n_handles * 3 > get_fileno_limit()) {
		int o_fnl = get_fileno_limit();
		increase_nofile_limit(kroundup32(settings.n_handles));
		fprintf(stderr, "[%s] Increased nofile limit from %i to %i.\n", __func__, o_fnl,
				kroundup32(settings.n_handles));
	}
	if(settings.ffq_prefix && !settings.run_hash_dmp) {
		fprintf(stderr, "[E:%s] Final fastq prefix option provided but run_hash_dmp not selected."
				"Either eliminate the -f flag or add the -d flag.\n", __func__);
		exit(EXIT_FAILURE);
	}

	if(!settings.homing_sequence) {
		fprintf(stderr, "[E:%s] Homing sequence not provided. Required.\n", __func__);
		exit(EXIT_FAILURE);
	}
	else {
#if !NDEBUG
		fprintf(stderr, "[D:%s] Homing sequence: %s.\n", __func__, settings.homing_sequence);
#endif
		for(int i = 0; settings.homing_sequence[i]; ++i) {
			switch(settings.homing_sequence[i]) {
			case 'A': break;
			case 'C': break;
			case 'G': break;
			case 'T': break;
			case 'a': // Fall-through
			case 'g': // Fall-through
			case 'c': // Fall-through
			case 't': settings.homing_sequence[i] -= 32; // Converts lower-case to upper-case
            default: fprintf(stderr, "[E:%s] Homing sequence contains illegal characters. Accepted: ACGT. Parameter: %s.\n",
                             __func__, settings.homing_sequence);
			exit(EXIT_FAILURE);
			}
		}
	}
	if(!settings.blen) {
		fprintf(stderr, "[E:%s] Barcode length not provided. Required. Abort!\n", __func__);
		exit(EXIT_FAILURE);
	}

	if(settings.max_blen > 0 && settings.max_blen * 2 < settings.blen) {
		fprintf(stderr, "[E:%s] Max blen (%i) must be less than the minimum blen provided (%i).\n",
                __func__, settings.max_blen, settings.blen / 2);
        exit(EXIT_FAILURE);
	}

	if(settings.offset)
		settings.blen -= 2 * settings.offset;


	if(argc - 1 != optind + 1) {
		fprintf(stderr, "[E:%s] Both read 1 and read 2 fastqs are required. See usage.\n", __func__);
		print_crms_usage(argv[0]);
		return 1;
	}
	settings.input_r1_path = strdup(argv[optind]);
	settings.input_r2_path = strdup(argv[optind + 1]);

	if(!settings.output_basename) {
		settings.output_basename = (char *)malloc(21 * sizeof(char));
		rand_string(settings.output_basename, 20);
		fprintf(stderr, "[%s] Output basename not provided. Defaulting to random: %s.\n",
                __func__, settings.output_basename);
	}
	mark_splitter_t *splitter = (settings.annealed) ? pp_split_annealed(&settings): pp_split_inline(&settings);
	cond_free(settings.rescaler);
	if(settings.run_hash_dmp) {
		fprintf(stderr, "[%s] Now executing hashmap-powered read collapsing and molecular demultiplexing.\n",
                __func__);
		if(!settings.ffq_prefix) {
            int has_period = 0;
            for(int i = 0; settings.input_r1_path[i]; ++i) {
                if(settings.input_r1_path[i] == '.') {
                    has_period = 1; break;
                }
            }
            if(has_period) {
			    settings.ffq_prefix = make_default_outfname(settings.input_r1_path, ".dmp.final");
                fprintf(stderr, "[%s] No output final prefix set. Defaulting to variation on input ('%s').\n",
                        __func__, settings.ffq_prefix);
            }
            else {
                settings.ffq_prefix = (char *)malloc(21 * sizeof(char));
                rand_string(settings.ffq_prefix, 20);
                fprintf(stderr, "[%s] No output final prefix set. Selecting random output name ('%s').\n",
                        __func__, settings.ffq_prefix);
            }
		}
		// Whatever I end up putting into here.
		splitterhash_params_t *params = init_splitterhash(&settings, splitter);
		char del_buf[500];
		#pragma omp parallel
		{
			#pragma omp for schedule(dynamic, 1)
			for(int i = 0; i < settings.n_handles; ++i) {
				fprintf(stderr, "[%s] Now running hash dmp core on input filename %s and output filename %s.\n",
						__func__, params->infnames_r1[i], params->outfnames_r1[i]);
				stranded_hash_dmp_core(params->infnames_r1[i], params->outfnames_r1[i]);
			}
		}
			#pragma omp for schedule(dynamic, 1)
			for(int i = 0; i < settings.n_handles; ++i) {
				fprintf(stderr, "[%s] Now running hash dmp core on input filename "
                        "%s and output filename %s.\n",
						__func__, params->infnames_r2[i], params->outfnames_r2[i]);
				stranded_hash_dmp_core(params->infnames_r2[i], params->outfnames_r2[i]);
			}

		if(settings.cleanup) {
	        #pragma omp parallel for schedule(dynamic, 1)
	         for(int i = 0; i < settings.n_handles; ++i) {
				char tmpbuf[500];
				fprintf(stderr, "[%s] Now removing temporary files '%s', '%s'.\n",
			            __func__, params->infnames_r1[i], params->infnames_r2[i]);
				sprintf(tmpbuf, "rm %s %s", params->infnames_r2[i], params->infnames_r1[i]);
	            CHECK_CALL(tmpbuf);
			}
        }
		// Remove temporary split files
        //
		char cat_buff[CAT_BUFFER_SIZE];
		char ffq_r1[200];
		char ffq_r2[200];
		sprintf(ffq_r1, "%s.R1.fq", settings.ffq_prefix);
		sprintf(ffq_r2, "%s.R2.fq", settings.ffq_prefix);
		sprintf(cat_buff, settings.gzip_output ? "> %s.gz" : "> %s", ffq_r1);
		CHECK_CALL(cat_buff);
		sprintf(cat_buff, settings.gzip_output ? "> %s.gz" : "> %s", ffq_r2);
		CHECK_CALL(cat_buff);
		char cat_buff1[CAT_BUFFER_SIZE] = "/bin/cat ";
		if(settings.panthera) {
			char cat_buff2[CAT_BUFFER_SIZE] = "/bin/cat ";
			for(int i = 0; i < settings.n_handles; ++i) {
				strcat(cat_buff1, params->outfnames_r1[i]);
				strcat(cat_buff1, " ");
				strcat(cat_buff2, params->outfnames_r2[i]);
				strcat(cat_buff2, " ");
			}
			if(settings.gzip_output) {
				sprintf(del_buf, " | pigz - -%i -p %i", settings.gzip_compression, settings.threads >= 2 ? settings.threads >> 1: 1);
				strcat(cat_buff1, del_buf);
				strcat(cat_buff2, del_buf);
			}
			strcat(cat_buff1, " > ");
			strcat(cat_buff1, ffq_r1);
			strcat(cat_buff2, " > ");
			strcat(cat_buff2, ffq_r2);
			if(settings.gzip_output) {
				strcat(cat_buff1, ".gz");
				strcat(cat_buff2, ".gz");
			}
			FILE *c1_popen = popen(cat_buff1, "w");
			FILE *c2_popen = popen(cat_buff2, "w");
			if(pclose(c2_popen) || pclose(c1_popen)) {
				fprintf(stderr, "[%s] Background cat command failed. ('%s').\n", __func__, cat_buff1);
				exit(EXIT_FAILURE);
			}
		}
		else {
			for(int i = 0; i < settings.n_handles; ++i) {
				// Clear files if present
				if(settings.gzip_output)
					sprintf(cat_buff, "cat %s | pigz - -p %i -%i >> %s.gz",
							params->outfnames_r1[i], settings.threads >= 2 ? settings.threads >> 1: 1, settings.gzip_compression, ffq_r1);
				else
					sprintf(cat_buff, "cat %s >> %s", params->outfnames_r1[i], ffq_r1);
				FILE *g1_popen = popen(cat_buff, "w");
				if(settings.gzip_output)
					sprintf(cat_buff, "cat %s | pigz - -p %i -%i >> %s.gz",
							params->outfnames_r2[i], settings.threads >= 2 ? settings.threads >> 1: 1, settings.gzip_compression, ffq_r2);
				else
					sprintf(cat_buff, "cat %s >> %s", params->outfnames_r2[i], ffq_r2);
				FILE *g2_popen = popen(cat_buff, "w");
				if(pclose(g2_popen) || pclose(g1_popen)){
					fprintf(stderr, "[E:%s] Background system call failed.\n", __func__);
					exit(EXIT_FAILURE);
				}
			}
		}
		if(settings.cleanup) {
			#pragma omp parallel for shared(params)
			for(int i = 0; i < params->n; ++i) {
				char tmpbuf[500];
				sprintf(tmpbuf, "rm %s %s", params->outfnames_r1[i], params->outfnames_r2[i]);
				CHECK_CALL(tmpbuf);
			}
		}
		splitterhash_destroy(params);
		free(settings.ffq_prefix);
	}
    free_marksplit_settings(settings);
    splitter_destroy(splitter);
	return 0;
}

inline int test_homing_seq(kseq_t *seq1, kseq_t *seq2, marksplit_settings_t *settings_ptr)
{
    return (settings_ptr->homing_sequence) ?
        (memcmp(seq1->seq.s + (settings_ptr->blen / 2 + settings_ptr->offset),
               settings_ptr->homing_sequence, settings_ptr->homing_sequence_length) == 0): 1;
}

inline char test_hp_inline(char *barcode, int length, int threshold)
{
	int run = 0;
	char last = 0;
	for(int i = 0; i < length; i++){
		if(barcode[i] == 'N')
			return '0';
		if(barcode[i] == last) {
			if(++run >= threshold)
                return '0';
		}
		else
			run = 0, last = barcode[i];
	}
	return '1';
}

