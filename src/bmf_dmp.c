/*
 * Conditional reverse complement Rescale Mark Split
 * For inline Loeb adapters only.
 */

#include "bmf_dmp.h"


void print_crms_usage(char *executable)
{
		fprintf(stderr, "Usage: %s <options> <Fq.R1.seq> <Fq.R2.seq>"
						"\nFlags:\n"
						"-$: Flag for single-end mode.\n"
						"-=: Emit interleaved final output to stdout.\n"
						"-l: Number of nucleotides at the beginning of each read to "
						"use for barcode. Final barcode length is twice this. REQUIRED.\n"
						"-s: homing sequence. REQUIRED.\n"
						"-o: Mark/split temporary file basename. Defaults to a random string.\n"
						"-t: Homopolymer failure threshold. A molecular barcode with"
						" a homopolymer of length >= this limit is flagged as QC fail."
						"Default: 10.\n"
						"-n: Number of nucleotides at the beginning of the barcode to use to split the output. Default: 4.\n"
						"-m: Mask first n nucleotides in read for barcode. Default: 0.\n"
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
	LOG_ERROR("Unrecognized option %s for flag %c. Abort!\n", optarg, optopt);
}

void make_outfname(marksplit_settings_t *settings)
{
	int has_period = 0;
	for(int i = 0; settings->input_r1_path[i]; ++i) {
		if(settings->input_r1_path[i] == '.') {
			has_period = 1; break;
		}
	}
	if(has_period) {
	settings->ffq_prefix = make_default_outfname(settings->input_r1_path, ".dmp.final");
		fprintf(stderr, "[%s] No output final prefix set. Defaulting to variation on input ('%s').\n",
				__func__, settings->ffq_prefix);
	} else {
		settings->ffq_prefix = (char *)malloc((RANDSTR_SIZE + 1) * sizeof(char));
		rand_string(settings->ffq_prefix, RANDSTR_SIZE);
		fprintf(stderr, "[%s] No output final prefix set. Selecting random output name ('%s').\n",
				__func__, settings->ffq_prefix);
	}
}

void cleanup_hashdmp(marksplit_settings_t *settings, splitterhash_params_t *params)
{
	if(settings->cleanup) {
		#pragma omp parallel for
		for(int i = 0; i < params->n; ++i) {
			char tmpbuf[1000];
			if(settings->is_se) sprintf(tmpbuf, "rm %s", params->outfnames_r1[i]);
			else sprintf(tmpbuf, "rm %s %s", params->outfnames_r1[i], params->outfnames_r2[i]);
			CHECK_CALL(tmpbuf);
		}
	}
}


void parallel_hash_dmp_core(marksplit_settings_t *settings, splitterhash_params_t *params, hash_dmp_fn func)
{
	#pragma omp parallel for schedule(dynamic, 1)
	for(int i = 0; i < settings->n_handles; ++i) {
		fprintf(stderr, "[%s] Now running hash dmp core on input filename %s and output filename %s.\n",
				__func__, params->infnames_r1[i], params->outfnames_r1[i]);
		func(params->infnames_r1[i], params->outfnames_r1[i], settings->gzip_compression);
		if(settings->cleanup) {
			char tmpbuf[500];
			sprintf(tmpbuf, "rm %s", params->infnames_r1[i]);
			CHECK_CALL(tmpbuf);
		}
	}
	if(settings->is_se) return;
	#pragma omp parallel for schedule(dynamic, 1)
	for(int i = 0; i < settings->n_handles; ++i) {
		fprintf(stderr, "[%s] Now running hash dmp core on input filename %s and output filename %s.\n",
				__func__, params->infnames_r2[i], params->outfnames_r2[i]);
		func(params->infnames_r2[i], params->outfnames_r2[i], settings->gzip_compression);
		kstring_t ks = {0, 0, NULL};
		if(settings->cleanup) {
			ksprintf(&ks, "rm %s", params->infnames_r2[i]);
			CHECK_CALL(ks.s);
			ks.l = 0;
		}
		free(ks.s);
	}
}


void call_clowder_se(marksplit_settings_t *settings, splitterhash_params_t *params, char *ffq_r1)
{
	// Clear output files.
	kstring_t ks = {0, 0, NULL};
	ksprintf(&ks, settings->gzip_output ? "> %s.gz" : "> %s", ffq_r1);
	CHECK_CALL(ks.s);
	for(int i = 0; i < settings->n_handles; ++i) {
		ks.l = 0;
		// Clear files if present
		ksprintf(&ks, "cat %s >> %s", params->outfnames_r1[i], ffq_r1);
		CHECK_POPEN(ks.s);
	}
	free(ks.s);
}

void call_clowder_pe(marksplit_settings_t *settings, splitterhash_params_t *params, char *ffq_r1, char *ffq_r2)
{
	// Clear output files.
	kstring_t ks = {0, 0, NULL};
	ksprintf(&ks, settings->gzip_output ? "> %s.gz" : "> %s", ffq_r1);
	CHECK_CALL(ks.s);
	ks.l = 0;
	ksprintf(&ks, settings->gzip_output ? "> %s.gz" : "> %s", ffq_r2);
	CHECK_CALL(ks.s);
	for(int i = 0; i < settings->n_handles; ++i) {
		// Clear files if present
		ks.l = 0;
		ksprintf(&ks, "cat %s >> %s", params->outfnames_r1[i], ffq_r1);
		FILE *g1_popen = popen(ks.s, "w");
		ks.l = 0;
		ksprintf(&ks, "cat %s >> %s", params->outfnames_r2[i], ffq_r2);
		FILE *g2_popen = popen(ks.s, "w");
		if(pclose(g2_popen) || pclose(g1_popen)){
			LOG_ERROR("Background system call failed.\n");
		}
	}
	free(ks.s);
}


void call_panthera_se(marksplit_settings_t *settings, splitterhash_params_t *params, char *ffq_r1)
{
	kstring_t ks = {0, 0, NULL};
	// Clear output files.
	ksprintf(&ks, settings->gzip_output ? "> %s.gz" : "> %s", ffq_r1);
	CHECK_CALL(ks.s);
	ks.l = 0;
	ksprintf(&ks, "/bin/cat ");
	for(int i = 0; i < settings->n_handles; ++i) {
		if(!isfile(params->outfnames_r1[i])) {
			LOG_ERROR("Output filename is not a file. Abort! ('%s').\n", params->outfnames_r1[i]);
		}
		ksprintf(&ks, " %s", params->outfnames_r1[i]);
	}
	ksprintf(&ks, " > %s", ffq_r1);
	if(settings->gzip_output) kputs(".gz", &ks);
	CHECK_POPEN(ks.s);
	free(ks.s);
}

void check_rescaler(marksplit_settings_t *settings, int arr_size)
{
	if(settings->rescaler) {
		for(int i = 0; i < arr_size; ++i) {
			if(settings->rescaler[i] < 0) {
				LOG_ERROR("Rescaler's got a negative number in pp_split_inline."
						" %i. Index: %i.\n", settings->rescaler[i], i);
			}
			else if(settings->rescaler[i] == 0) {
				LOG_ERROR("Rescaler's got a zero in pp_split_inline."
						" %i. Index: %i.\n", settings->rescaler[i], i);
			}
		}
	}
}

void call_stdout(marksplit_settings_t *settings, splitterhash_params_t *params, char *ffq_r1, char *ffq_r2)
{
	char fname_buf[500];
	kstring_t str1 = {0, 0, NULL}, str2 = {0, 0, NULL};
	kputs("cat ", &str1);
	ks_resize(&str1, 1 << 16);
	for(int i = 0; i < settings->n_handles; ++i) {
		sprintf(fname_buf, " %s", params->outfnames_r1[i]);
		while(str1.m < strlen(fname_buf) + str1.l + 1) ks_resize(&str1, str1.m << 1);
		strcat(str1.s, fname_buf); str1.l += strlen(fname_buf);
	}
	const char suffix[] = " | paste -d'~' - - - - ";
	while(str1.m < sizeof(suffix) + str1.l) ks_resize(&str1, str1.m << 1);
	strcat(str1.s, suffix);
	str1.l += sizeof(suffix);
	str2 = str1; // Copy everything, including a pointer that str2 doesn't own.
	str2.s = strdup(str1.s); // strdup the string.
	for(uint32_t i = 0; i < str2.l - 1; ++i)
		if(str2.s[i] == 'R' && str2.s[i + 1] == '1')
			str2.s[i + 1] = '2';

	const char final_template[] = "pr -mts'~' <(%s) <(%s) | tr '~' '\\n'";
	char *final = (char *)malloc(str1.m + str2.m + sizeof(final_template) / sizeof(char)); // Should be plenty of space
	sprintf(final, final_template, str1.s, str2.s);
	bash_system(final);
	free(str1.s), free(str2.s);
	free(final);
}


void call_clowder(marksplit_settings_t *settings, splitterhash_params_t *params, char *ffq_r1, char *ffq_r2)
{
	settings->is_se ? call_clowder_se(settings, params, ffq_r1):
			call_clowder_pe(settings, params, ffq_r1, ffq_r2);
}

void call_panthera(marksplit_settings_t *settings, splitterhash_params_t *params, char *ffq_r1, char *ffq_r2)
{
	settings->is_se ? call_panthera_se(settings, params, ffq_r1):
			call_panthera_pe(settings, params, ffq_r1, ffq_r2);
}

void call_panthera_pe(marksplit_settings_t *settings, splitterhash_params_t *params, char *ffq_r1, char *ffq_r2)
{
	kstring_t ks1 = {0, 0, NULL};
	kputs("> ", &ks1), kputs(ffq_r1, &ks1);
	if(settings->gzip_output) kputs(".gz", &ks1);
	// Clear output files.
	//ksprintf(&ks1, settings->gzip_output ? "> %s.gz" : "> %s", ffq_r1);
	CHECK_CALL(ks1.s); ks1.l = 0;
	ksprintf(&ks1, settings->gzip_output ? "> %s.gz" : "> %s", ffq_r2);
	CHECK_CALL(ks1.s); ks1.l = 0;
	kputs("/bin/cat ", &ks1);
	kstring_t ks2 = {0};
	ksprintf(&ks2, ks1.s);
	for(int i = 0; i < settings->n_handles; ++i) {
		if(!isfile(params->outfnames_r1[i])) {
			LOG_ERROR("Output filename is not a file. Abort! ('%s').\n", params->outfnames_r1[i]);
		}
		if(!isfile(params->outfnames_r2[i])) {
			LOG_ERROR("Output filename is not a file. Abort! ('%s').\n", params->outfnames_r2[i]);
		}
		ksprintf(&ks1, "%s ", params->outfnames_r1[i]);
		ksprintf(&ks2, "%s ", params->outfnames_r2[i]);
	}
	ksprintf(&ks1, " > %s", ffq_r1);
	ksprintf(&ks2, " > %s", ffq_r2);
	if(settings->gzip_output)
		kputs(".gz", &ks1), kputs(".gz", &ks2);
	FILE *c1_popen = popen(ks1.s, "w");
	FILE *c2_popen = popen(ks2.s, "w");
	if(pclose(c2_popen) || pclose(c1_popen)) {
		LOG_ERROR("Background cat command failed. ('%s' or '%s').\n", ks1.s, ks2.s);
	}
	free(ks1.s), free(ks2.s);
}

void clean_homing_sequence(char *sequence) {
	while(*sequence) {
		switch(*sequence) {
		case 'A': break;
		case 'C': break;
		case 'G': break;
		case 'T': break;
		case 'a': // Fall-through
		case 'g': // Fall-through
		case 'c': // Fall-through
		case 't': *sequence -= UPPER_LOWER_OFFSET; break;// Converts lower-case to upper-case
		default:
			LOG_ERROR("Homing sequence contains illegal characters. Accepted: [acgtACGT]. Character: %c.\n",
					*sequence);
		}
		++sequence;
	}
}

/*
 * Pre-processes (pp) and splits fastqs with inline barcodes.
mark_splitter_t *pp_split_inline_se(marksplit_settings_t *settings)
{
	fprintf(stderr, "[%s] Opening fastq file '%s'.\n", __func__, settings->input_r1_path);
	if(!isfile(settings->input_r1_path)) {
		LOG_ERROR("Could not open read paths: at least one is not a file.\n");
	}

	if(settings->rescaler_path) settings->rescaler = parse_1d_rescaler(settings->rescaler_path);

	mark_splitter_t *splitter = (mark_splitter_t *)malloc(sizeof(mark_splitter_t));
	*splitter = init_splitter(settings);
	gzFile fp = gzopen(settings->input_r1_path, "r");
	kseq_t *seq = kseq_init(fp);
	int l;
	size_t count = 0;
	int pass_fail = 1;
	if((l = kseq_read(seq)) < 0) {
		free_marksplit_settings(*settings);
		splitter_destroy(splitter);
		LOG_ERROR("Could not open fastq for reading. Abort!\n");
	}
	LOG_DEBUG("Read length (inferred): %lu.\n", seq->seq.l);
	check_rescaler(settings, seq->seq.l * 4 * 2 * nqscores);
	tmp_mseq_t *tmp = init_tm_ptr(seq->seq.l, settings->blen);
	const int default_nlen = settings->blen + settings->offset + settings->homing_sequence_length;
	int n_len = nlen_homing_se(seq, settings, default_nlen, &pass_fail);
	mseq_t *rseq = mseq_rescale_init(seq, settings->rescaler, tmp, 0);
	rseq->barcode[settings->blen] = '\0';
	pass_fail &= test_hp(rseq->barcode, settings->hp_threshold); // Fail if test_hp is 0.
	memcpy(rseq->barcode, seq->seq.s + settings->offset, settings->blen);
	// Get first barcode.
	update_mseq(rseq, seq, settings->rescaler, tmp, n_len, 0);
	uint64_t bin = get_binner_type(rseq->barcode, settings->n_nucs, uint64_t);
	mseq2fq(splitter->tmp_out_handles_r1[bin], rseq, pass_fail, rseq->barcode);
	while (LIKELY((l = kseq_read(seq)) >= 0)) {
		if(UNLIKELY(++count % settings->notification_interval == 0)) {
			LOG_INFO("Number of records processed: %lu.\n", count);
		}
		// Iterate through second fastq file.
		n_len = nlen_homing_se(seq, settings, default_nlen, &pass_fail);
		update_mseq(rseq, seq, settings->rescaler, tmp, n_len, 0);
		memcpy(rseq->barcode, seq->seq.s + settings->offset, settings->blen);
		pass_fail &= test_hp(rseq->barcode, settings->hp_threshold);
		mseq2fq(splitter->tmp_out_handles_r1[bin], rseq, pass_fail, rseq->barcode);
	}
	for(int i = 0; i < splitter->n_handles; ++i) fclose(splitter->tmp_out_handles_r1[i]);
	tm_destroy(tmp);
	mseq_destroy(rseq);
	kseq_destroy(seq);
	gzclose(fp);
	return splitter;
}
 */

/*
 * Pre-processes (pp) and splits fastqs with inline barcodes.
 */
mark_splitter_t *pp_split_inline(marksplit_settings_t *settings)
{
#if WRITE_BARCODE_FQ
	gzFile bcfp1 = gzopen("tmp.molbc.r1.fq.gz", "wb4");
	gzFile bcfp2 = gzopen("tmp.molbc.r2.fq.gz", "wb4");
#endif
	LOG_INFO("Opening fastq files %s and %s.\n", settings->input_r1_path, settings->input_r2_path);
	if(!(strcmp(settings->input_r1_path, settings->input_r2_path))) {
		LOG_ERROR("Hey, it looks like you're trying to use the same path for both r1 and r2. "
				"At least try to fool me by making a symbolic link.\n");
	}
	if(!isfile(settings->input_r1_path) || !isfile(settings->input_r2_path)) {
		LOG_ERROR("Could not open read paths: at least one is not a file.\n");
	}
	if(settings->rescaler_path) settings->rescaler = parse_1d_rescaler(settings->rescaler_path);
	mark_splitter_t *splitter = (mark_splitter_t *)malloc(sizeof(mark_splitter_t));
	*splitter = init_splitter(settings);
	gzFile fp1 = gzopen(settings->input_r1_path, "r");
	gzFile fp2 = gzopen(settings->input_r2_path, "r");
	kseq_t *seq1 = kseq_init(fp1), *seq2 = kseq_init(fp2);
	int l1, l2;
	int pass_fail = 1;
	if((l1 = kseq_read(seq1)) < 0 || (l2 = kseq_read(seq2)) < 0) {
			free_marksplit_settings(*settings);
			splitter_destroy(splitter);
			LOG_ERROR("Could not open fastqs for reading. Abort!\n");
	}
	LOG_DEBUG("Read length (inferred): %lu.\n", seq1->seq.l);
	check_rescaler(settings, seq1->seq.l * 4 * 2 * nqscores);
	tmp_mseq_t *tmp = init_tm_ptr(seq1->seq.l, settings->blen);
	int switch_reads = switch_test(seq1, seq2, settings->offset);
	const int default_nlen = settings->blen1_2 + settings->offset + settings->homing_sequence_length;
	int n_len = nlen_homing_default(seq1, seq2, settings, default_nlen, &pass_fail);
	mseq_t *rseq1, *rseq2;
	rseq1 = mseq_rescale_init(seq1, settings->rescaler, tmp, 0);
	rseq2 = mseq_rescale_init(seq2, settings->rescaler, tmp, 1);
	rseq1->barcode[settings->blen] = '\0';
#if WRITE_BARCODE_FQ
	write_bc_to_file(bcfp1, bcfp2, seq1, seq2, settings);
#endif
	if(switch_reads) {
		memcpy(rseq1->barcode, seq2->seq.s + settings->offset, settings->blen1_2);
		memcpy(rseq1->barcode + settings->blen1_2, seq1->seq.s + settings->offset, settings->blen1_2);
	} else {
		memcpy(rseq1->barcode, seq1->seq.s + settings->offset, settings->blen1_2);
		memcpy(rseq1->barcode + settings->blen1_2, seq2->seq.s + settings->offset, settings->blen1_2);
	}
	pass_fail &= test_hp(rseq1->barcode, settings->hp_threshold);
	mask_mseq(rseq1, n_len); mask_mseq(rseq2, n_len);
	// Get first barcode.
	update_mseq(rseq1, seq1, settings->rescaler, tmp, n_len, 0);
	update_mseq(rseq2, seq2, settings->rescaler, tmp, n_len, 1);
	uint64_t bin = get_binner_type(rseq1->barcode, settings->n_nucs, uint64_t);
	assert(bin < settings->n_handles);
	if(switch_reads) {
		mseq2fq_stranded(splitter->tmp_out_handles_r1[bin], rseq2, pass_fail, rseq1->barcode, 'R');
		mseq2fq_stranded(splitter->tmp_out_handles_r2[bin], rseq1, pass_fail, rseq1->barcode, 'R');
	} else {
		mseq2fq_stranded(splitter->tmp_out_handles_r1[bin], rseq1, pass_fail, rseq1->barcode, 'F');
		mseq2fq_stranded(splitter->tmp_out_handles_r2[bin], rseq2, pass_fail, rseq1->barcode, 'F');
	}
	uint64_t count = 0;
	while (LIKELY((l1 = kseq_read(seq1)) >= 0) && LIKELY((l2 = kseq_read(seq2)) >= 0)) {
		if(UNLIKELY(++count % settings->notification_interval == 0))
			LOG_INFO("Number of records processed: %lu.\n", count);
		// Sets pass_fail
		n_len = nlen_homing_default(seq1, seq2, settings, default_nlen, &pass_fail);
#if WRITE_BARCODE_FQ
		write_bc_to_file(bcfp1, bcfp2, seq1, seq2, settings);
#endif
		// Update mseqs
		update_mseq(rseq1, seq1, settings->rescaler, tmp, n_len, 0);
		update_mseq(rseq2, seq2, settings->rescaler, tmp, n_len, 1);

		if(switch_test(seq1, seq2, settings->offset)) {
			// Copy barcode over
			memcpy(rseq1->barcode, seq2->seq.s + settings->offset, settings->blen1_2);
			memcpy(rseq1->barcode + settings->blen1_2, seq1->seq.s + settings->offset, settings->blen1_2);
			// Test for homopolymer failure
			pass_fail &= test_hp(rseq1->barcode, settings->hp_threshold);
			bin = get_binner_type(rseq1->barcode, settings->n_nucs, uint64_t);
			assert(bin < settings->n_handles);
			// Write out
			mseq2fq_stranded(splitter->tmp_out_handles_r1[bin], rseq2, pass_fail, rseq1->barcode, 'R');
			mseq2fq_stranded(splitter->tmp_out_handles_r2[bin], rseq1, pass_fail, rseq1->barcode, 'R');
		} else {
			// Copy barcode over
			memcpy(rseq1->barcode, seq1->seq.s + settings->offset, settings->blen1_2);
			memcpy(rseq1->barcode + settings->blen1_2, seq2->seq.s + settings->offset, settings->blen1_2);
			// Test for homopolymer failure
			pass_fail &= test_hp(rseq1->barcode, settings->hp_threshold);
			bin = get_binner_type(rseq1->barcode, settings->n_nucs, uint64_t);
			assert(bin < (uint64_t)settings->n_handles);
			// Write out
			mseq2fq_stranded(splitter->tmp_out_handles_r1[bin], rseq1, pass_fail, rseq1->barcode, 'F');
			mseq2fq_stranded(splitter->tmp_out_handles_r2[bin], rseq2, pass_fail, rseq1->barcode, 'F');
		}
	}
	LOG_DEBUG("Cleaning up.\n");
	for(l1 = 0; l1 < splitter->n_handles; ++l1)
		fclose(splitter->tmp_out_handles_r1[l1]), fclose(splitter->tmp_out_handles_r2[l1]);
	tm_destroy(tmp);
	mseq_destroy(rseq1), mseq_destroy(rseq2);
	kseq_destroy(seq1), kseq_destroy(seq2);
	gzclose(fp1), gzclose(fp2);
#ifdef WRITE_BARCODE_FQ
	gzclose(bcfp1), gzclose(bcfp2);
#endif
	return splitter;
}


int dmp_main(int argc, char *argv[])
{
	if(argc == 1) print_crms_usage(argv[0]), exit(EXIT_FAILURE);
	// Build settings struct

	marksplit_settings_t settings = {0};

	settings.hp_threshold = 10;
	settings.n_nucs = 2;
	settings.notification_interval = 1000000;
	settings.threads = 1;
	settings.max_blen = -1;
	settings.gzip_compression = 1;
	settings.cleanup = 1;

	//omp_set_dynamic(0); // Tell omp that I want to set my number of threads 4realz
	int c;
	while ((c = getopt(argc, argv, "t:o:n:s:l:m:r:p:f:v:u:g:i:zwcdh?$=")) > -1) {
		switch(c) {
			case 'c': settings.panthera = 1; break;
			case 'd': settings.run_hash_dmp = 1; break;
			case 'f': settings.ffq_prefix = strdup(optarg); break;
			case 'g': settings.gzip_compression = atoi(optarg); break;
			case 'l': settings.blen = atoi(optarg); break;
			case 'm': settings.offset = atoi(optarg); break;
			case 'n': settings.n_nucs = atoi(optarg); break;
			case 'o': settings.tmp_basename = strdup(optarg); break;
			case 'p': settings.threads = atoi(optarg); break;
			case 'r': settings.rescaler_path = strdup(optarg); break;
			case 's': settings.homing_sequence = strdup(optarg); settings.homing_sequence_length = strlen(settings.homing_sequence); break;
			case 't': settings.hp_threshold = atoi(optarg); break;
			case 'u': settings.notification_interval = atoi(optarg); break;
			case 'v': settings.max_blen = atoi(optarg); break;
			case 'w': settings.cleanup = 0; break;
			case 'z': settings.gzip_output = 1; break;
			case '$': settings.is_se = 1; break;
			case '=': settings.to_stdout = 1; break;
			case '?': // Fall-through
			case 'h': print_crms_usage(argv[0]), exit(EXIT_SUCCESS);
		}
	}

	if(!settings.gzip_output) settings.gzip_compression = 0;

	// Check for proper command-line usage.
	if(settings.is_se) {
		if(argc < 4) print_crms_usage(argv[0]), exit(EXIT_FAILURE);
		if(argc != optind + 1) {
			fprintf(stderr, "[E:%s] Exactly one read fastq required for single-end. See usage.\n", __func__);
			print_crms_usage(argv[0]);
			return EXIT_FAILURE;
		}
		// Number of file handles
		settings.n_handles = ipow(4, settings.n_nucs);
		if(settings.n_handles * 2 > get_fileno_limit()) {
			increase_nofile_limit(settings.n_handles * 2);
			fprintf(stderr, "[%s] Increased nofile limit from %i to %i.\n", __func__, get_fileno_limit(),
					settings.n_handles * 2);
		}
		// Handle filenames
		settings.input_r1_path = strdup(argv[optind]);
	} else {
		if(argc < 5) print_crms_usage(argv[0]), exit(EXIT_FAILURE);
		if(argc != optind + 2) {
			fprintf(stderr, "[E:%s] Both read 1 and read 2 fastqs are required for paired-end. See usage.\n", __func__);
			print_crms_usage(argv[0]);
			return EXIT_FAILURE;
		}
		// Number of file handles
		settings.n_handles = ipow(4, settings.n_nucs);
		if(settings.n_handles * 4 > get_fileno_limit()) {
			increase_nofile_limit(settings.n_handles * 4);
			fprintf(stderr, "[%s] Increased nofile limit from %i to %i.\n", __func__, get_fileno_limit(),
					settings.n_handles * 4);
		}
		// Handle filenames
		settings.input_r1_path = strdup(argv[optind]);
		settings.input_r2_path = strdup(argv[optind + 1]);
	}
	// Required parameters
	if(settings.ffq_prefix && !settings.run_hash_dmp) {
		fprintf(stderr, "[E:%s] Final fastq prefix option provided but run_hash_dmp not selected."
				"Either eliminate the -f flag or add the -d flag.\n", __func__);
		exit(EXIT_FAILURE);
	}

	// Handle number of threads
	omp_set_num_threads(settings.threads);

	// Handle homing sequence
	if(!settings.homing_sequence) {
		fprintf(stderr, "[E:%s] Homing sequence not provided. Required.\n", __func__);
		exit(EXIT_FAILURE);
	}
#if !NDEBUG
	fprintf(stderr, "[D:%s] Homing sequence: %s.\n", __func__, settings.homing_sequence);
#endif
	clean_homing_sequence(settings.homing_sequence);

	// Handle barcode length
	if(!settings.blen) {
		fprintf(stderr, "[E:%s] Barcode length not provided. Required. Abort!\n", __func__);
		exit(EXIT_FAILURE);
	}
	// Handle the offset parameter. If false, blen doesn't change.
	if(settings.is_se) {
		settings.blen -= settings.offset;
		if(settings.max_blen > 0 && settings.max_blen < settings.blen) {
			fprintf(stderr, "[E:%s] Max blen (%i) must be less than the minimum blen provided (%i).\n",
					__func__, settings.max_blen, settings.blen / 2);
			exit(EXIT_FAILURE);
		}
		if(settings.max_blen < 0) settings.max_blen = settings.blen;
	} else {
		settings.blen = (settings.blen - settings.offset) * 2;
		if(settings.max_blen > 0 && settings.max_blen * 2 < settings.blen) {
			fprintf(stderr, "[E:%s] Max blen (%i) must be less than the minimum blen provided (%i).\n",
					__func__, settings.max_blen, settings.blen / 2);
			exit(EXIT_FAILURE);
		}
		settings.blen1_2 = settings.blen / 2;
		if(settings.max_blen < 0) settings.max_blen = settings.blen1_2;
	}

	if(!settings.tmp_basename) {
		// If tmp_basename unset, create a random temporary file prefix.
		settings.tmp_basename = (char *)malloc(21 * sizeof(char));
		rand_string(settings.tmp_basename, 20);
		fprintf(stderr, "[%s] Temporary basename not provided. Defaulting to random: %s.\n",
				__func__, settings.tmp_basename);
	}

	// Run core
	mark_splitter_t *splitter = pp_split_inline(&settings);
	//mark_splitter_t *splitter = settings.is_se ? pp_split_inline_se(&settings): pp_split_inline(&settings);
	splitterhash_params_t *params;
	char ffq_r1[500];
	char ffq_r2[500];
	if(!settings.run_hash_dmp) {
		fprintf(stderr, "[%s] mark/split complete.\n", __func__);
		goto cleanup;
	}
#if !NDEBUG
	fprintf(stderr, "[D:%s] Now executing hashmap-powered read collapsing and molecular demultiplexing.\n",
				__func__);
#endif
	if(!settings.ffq_prefix) make_outfname(&settings);
	params = init_splitterhash(&settings, splitter);
	// Run cores.
	parallel_hash_dmp_core(&settings, params, &stranded_hash_dmp_core);

	// Remove temporary split files.
	sprintf(ffq_r1, "%s.R1.fq", settings.ffq_prefix);
	sprintf(ffq_r2, "%s.R2.fq", settings.ffq_prefix);
	// Cat temporary files together.
	if(settings.to_stdout)
		call_stdout(&settings, params, ffq_r1, ffq_r2);
	else if(settings.panthera)
		call_panthera(&settings, params, ffq_r1, ffq_r2);
	else
		call_clowder(&settings, params, ffq_r1, ffq_r2);
	cleanup_hashdmp(&settings, params);
	splitterhash_destroy(params);
	cleanup:
	free_marksplit_settings(settings);
	splitter_destroy(splitter);
	fprintf(stderr, "[%s] Successfully completed bmftools dmp!\n", __func__);
	return EXIT_SUCCESS;
}

