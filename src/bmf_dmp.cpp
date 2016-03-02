/*
 * Conditional reverse complement Rescale Mark Split
 * For inline Loeb adapters only.
 */

#include "bmf_dmp.h"

kstring_t salted_rand_string(char *infname, size_t n_rand) {
    kstring_t ret = {0, 0, NULL};
    ksprintf(&ret, infname);
    char *tmp;
    /* Try to find the last of the string so that we salt the returned string with the input filename if there's a period.
     *
     */
    if((tmp = strrchr(ret.s, '.')) != NULL){
        *tmp = '\0';
    }
    ret.l = strlen(ret.s);
    ks_resize(&ret, ret.l + n_rand + 1);
    kputc('.', &ret);
    rand_string(ret.s + ret.l, n_rand);
    return ret;
}


void print_crms_usage(char *executable)
{
        fprintf(stderr, "Usage: bmftools %s <options> <Fq.R1.seq> <Fq.R2.seq>"
                        "\nFlags:\n"
                        "-S: Flag for single-end mode.\n"
                        "-=: Emit interleaved final output to stdout.\n"
                        "-l: Number of nucleotides at the beginning of each read to "
                        "use for barcode. Final barcode length is twice this. REQUIRED.\n"
                        "-s: homing sequence. REQUIRED.\n"
                        "-o: Mark/split temporary file basename. Defaults to a random string.\n"
                        "-t: Homopolymer failure threshold. A molecular barcode with"
                        " a homopolymer of length >= this limit is flagged as QC fail."
                        "Default: 10.\n"
                        "-n: Number of nucleotides at the beginning of the barcode to use to split the output. Default: %i.\n"
                        "-m: Mask first n nucleotides in read for barcode. Default: 0.\n"
                        "-p: Number of threads to use if running uthash_dmp.\n"
                        "-d: Use this flag to to run hash_dmp.\n"
                        "-f: If running hash_dmp, this sets the Final Fastq Prefix. \n"
                        "The Final Fastq files will be named '<ffq_prefix>.R1.fq' and '<ffq_prefix>.R2.fq'.\n"
                        "-r: Path to flat text file with rescaled quality scores. If not provided, it will not be used.\n"
                        "-v: Maximum barcode length for a variable length barcode dataset. If left as default value,"
                        " (-1), other barcode lengths will not be considered.\n"
                        "-z: Flag to write out final output as compressed. Default: False.\n"
                        "-T: If unset, write uncompressed plain text temporary files. If not, use that compression level for temporary files.\n"
                        "-g: Gzip compression ratio if writing gzipped. Default (if writing compressed): 1 (mostly to reduce I/O).\n"
                        "In addition, won't work for enormous filenames or too many arguments. Default: False.\n"
                        "-u: Set notification/update interval for split. Default: 1000000.\n"
                        "-w: Set flag to leave temporary files. Primarily for debugging.\n"
                        "-h: Print usage.\n",
					executable, DEFAULT_N_NUCS);

}

void print_crms_opt_err(char *arg, char *optarg, char optopt)
{
    print_crms_usage(arg);
    LOG_EXIT("Unrecognized option %s for flag %c. Abort!\n", optarg, optopt);
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
        kstring_t rs = salted_rand_string(settings->index_fq_path, RANDSTR_SIZE);
        settings->ffq_prefix = ks_release(&rs);
        LOG_INFO("No output final prefix set. Defaulting to variation on input ('%s').\n", settings->ffq_prefix);
    } else {
        char *tmp = (char *)malloc(RANDSTR_SIZE + 1);
        settings->ffq_prefix = rand_string(tmp, RANDSTR_SIZE);
    }
}

void cleanup_hashdmp(marksplit_settings_t *settings, splitterhash_params_t *params)
{
    if(!settings->cleanup) return;
    #pragma omp parallel for
    for(int i = 0; i < params->n; ++i) {
        kstring_t ks = {0, 0, NULL};
        ksprintf(&ks, "rm %s %s", params->outfnames_r1[i], settings->is_se ? "": params->outfnames_r2[i]);
        CHECK_CALL(ks.s);
        free(ks.s);
    }
}


void parallel_hash_dmp_core(marksplit_settings_t *settings, splitterhash_params_t *params, hash_dmp_fn func)
{
    #pragma omp parallel for schedule(dynamic, 1)
    for(int i = 0; i < settings->n_handles; ++i) {
        LOG_INFO("Now running hash dmp core on input filename %s and output filename %s.\n",
                params->infnames_r1[i], params->outfnames_r1[i]);
        func(params->infnames_r1[i], params->outfnames_r1[i], settings->gzip_compression);
        if(settings->cleanup) {
            kstring_t ks = {0, 0, NULL};
            ksprintf(&ks, "rm %s", params->infnames_r1[i]);
            CHECK_CALL(ks.s);
            free(ks.s);
        }
    }
    if(settings->is_se) return;
    #pragma omp parallel for schedule(dynamic, 1)
    for(int i = 0; i < settings->n_handles; ++i) {
        LOG_INFO("Now running hash dmp core on input filename %s and output filename %s.\n",
                params->infnames_r2[i], params->outfnames_r2[i]);
        func(params->infnames_r2[i], params->outfnames_r2[i], settings->gzip_compression);
        if(settings->cleanup) {
            kstring_t ks = {0, 0, NULL};
            ksprintf(&ks, "rm %s", params->infnames_r2[i]);
            CHECK_CALL(ks.s);
            free(ks.s);
        }
    }
}


void cat_fastqs_se(marksplit_settings_t *settings, splitterhash_params_t *params, char *ffq_r1)
{
    kstring_t ks = {0, 0, NULL};
    // Clear output files.
    ksprintf(&ks, settings->gzip_output ? "> %s.gz" : "> %s", ffq_r1);
    CHECK_CALL(ks.s);
    ks.l = 0;
    ksprintf(&ks, "/bin/cat ");
    for(int i = 0; i < settings->n_handles; ++i) {
        if(!isfile(params->outfnames_r1[i]))
            LOG_EXIT("Output filename is not a file. Abort! ('%s').\n", params->outfnames_r1[i]);
        ksprintf(&ks, " %s", params->outfnames_r1[i]);
    }
    ksprintf(&ks, " > %s", ffq_r1);
    if(settings->gzip_output) kputs(".gz", &ks);
    CHECK_POPEN(ks.s);
    free(ks.s);
}
/*
 * Make sure that no rescaler values are invalid
 */
void check_rescaler(marksplit_settings_t *settings, int arr_size)
{
    if(settings->rescaler)
        for(int i = 0; i < arr_size; ++i)
            if(settings->rescaler[i] <= 0)
                LOG_EXIT("Invalid value in rescaler %i at index %i.\n", settings->rescaler[i], i);
}
/*
 * Emits final results to stdout
 */

void call_stdout(marksplit_settings_t *settings, splitterhash_params_t *params, char *ffq_r1, char *ffq_r2)
{
    kstring_t str1 = {0, 0, NULL}, str2 = {0, 0, NULL};
    kstring_t final = {0, 0, NULL};
    kputs((settings->gzip_output)? "zcat": "cat", &str1);
    ks_resize(&str1, 1 << 16);
    for(int i = 0; i < settings->n_handles; ++i)
        ksprintf(&str1, " %s", params->outfnames_r1[i]);
    kputs(" | paste -d'~' - - - - ", &str1);
    str2.s = kstrdup(&str1); // strdup the string.
    for(uint32_t i = 0; i < str2.l; ++i) {
        LOG_DEBUG("Current str.s + i: %s.\n", str2.s + i);
        if(memcmp(str2.s + i, "R1", 2) == 0)
            str2.s[i + 1] = '2';
    }

    const char final_template[] = "pr -mts'~' <(%s) <(%s) | tr '~' '\\n'";
    ksprintf(&final, final_template, str1.s, str2.s);
    bash_system(final.s);
    free(str1.s), free(str2.s);
    free(final.s);
}

void cat_fastqs(marksplit_settings_t *settings, splitterhash_params_t *params, char *ffq_r1, char *ffq_r2)
{
    settings->is_se ? cat_fastqs_se(settings, params, ffq_r1):
            cat_fastqs_pe(settings, params, ffq_r1, ffq_r2);
}

void cat_fastqs_pe(marksplit_settings_t *settings, splitterhash_params_t *params, char *ffq_r1, char *ffq_r2)
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
            LOG_EXIT("Output filename is not a file. Abort! ('%s').\n", params->outfnames_r1[i]);
        }
        if(!isfile(params->outfnames_r2[i])) {
            LOG_EXIT("Output filename is not a file. Abort! ('%s').\n", params->outfnames_r2[i]);
        }
        ksprintf(&ks1, "%s ", params->outfnames_r1[i]);
        ksprintf(&ks2, "%s ", params->outfnames_r2[i]);
    }
    ksprintf(&ks1, " > %s", ffq_r1);
    ksprintf(&ks2, " > %s", ffq_r2);
    if(settings->gzip_output) {
        kputs(".gz", &ks1);
        kputs(".gz", &ks2);
    }
    FILE *c1_popen = popen(ks1.s, "w");
    FILE *c2_popen = popen(ks2.s, "w");
    if(pclose(c2_popen) || pclose(c1_popen)) {
        LOG_EXIT("Background cat command failed. ('%s' or '%s').\n", ks1.s, ks2.s);
    }
    free(ks1.s), free(ks2.s);
}

/*
 * Check for invalid characters and convert all lower-case to upper case.
 */
void clean_homing_sequence(char *sequence) {
    while(*sequence) {
        switch(*sequence) {
        case 'A': case 'C':  case 'G':  case 'T':
            break;
        case 'a': case 'g': case 'c': case 't':
            *sequence -= UPPER_LOWER_OFFSET; break;// Converts lower-case to upper-case
        default:
            LOG_EXIT("Homing sequence contains illegal characters. Accepted: [acgtACGT]. Character: %c.\n",
                     *sequence);
        }
        ++sequence;
    }
}

/*
 * Pre-processes (pp) and splits fastqs with inline barcodes.
 */
mark_splitter_t *pp_split_inline(marksplit_settings_t *settings)
{
    LOG_INFO("Opening fastq files %s and %s.\n", settings->input_r1_path, settings->input_r2_path);
    if(!(strcmp(settings->input_r1_path, settings->input_r2_path))) {
        LOG_EXIT("Hey, it looks like you're trying to use the same path for both r1 and r2. "
                "At least try to fool me by making a symbolic link.\n");
    }
    if(!isfile(settings->input_r1_path) || !isfile(settings->input_r2_path)) {
        LOG_EXIT("Could not open read paths: at least one is not a file.\n");
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
            LOG_EXIT("Could not open fastqs for reading. Abort!\n");
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
    if(switch_reads) {
        memcpy(rseq1->barcode, seq2->seq.s + settings->offset, settings->blen1_2);
        memcpy(rseq1->barcode + settings->blen1_2, seq1->seq.s + settings->offset, settings->blen1_2);
    } else {
        memcpy(rseq1->barcode, seq1->seq.s + settings->offset, settings->blen1_2);
        memcpy(rseq1->barcode + settings->blen1_2, seq2->seq.s + settings->offset, settings->blen1_2);
    }
    pass_fail &= test_hp(rseq1->barcode, settings->hp_threshold);
    // Get first barcode.
    update_mseq(rseq1, seq1, settings->rescaler, tmp, n_len, 0);
    update_mseq(rseq2, seq2, settings->rescaler, tmp, n_len, 1);
    uint64_t bin = get_binner_type(rseq1->barcode, settings->n_nucs, uint64_t);
    assert(bin < (uint64_t)settings->n_handles);
    if(switch_reads) {
        mseq2fq_stranded(splitter->tmp_out_handles_r1[bin], rseq2, pass_fail, rseq1->barcode, 'R');
        mseq2fq_stranded(splitter->tmp_out_handles_r2[bin], rseq1, pass_fail, rseq1->barcode, 'R');
    } else {
        mseq2fq_stranded(splitter->tmp_out_handles_r1[bin], rseq1, pass_fail, rseq1->barcode, 'F');
        mseq2fq_stranded(splitter->tmp_out_handles_r2[bin], rseq2, pass_fail, rseq1->barcode, 'F');
    }
    uint64_t count = 0;
    while (LIKELY(((l1 = kseq_read(seq1)) >= 0) && ((l2 = kseq_read(seq2)) >= 0))) {
        if(UNLIKELY(++count % settings->notification_interval == 0))
            LOG_INFO("Number of records processed: %lu.\n", count);
        // Sets pass_fail
        n_len = nlen_homing_default(seq1, seq2, settings, default_nlen, &pass_fail);
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
            assert(bin < (uint64_t)settings->n_handles);
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
        gzclose(splitter->tmp_out_handles_r1[l1]), gzclose(splitter->tmp_out_handles_r2[l1]);
    tm_destroy(tmp);
    mseq_destroy(rseq1), mseq_destroy(rseq2);
    kseq_destroy(seq1), kseq_destroy(seq2);
    gzclose(fp1), gzclose(fp2);
    return splitter;
}


int dmp_main(int argc, char *argv[])
{
    if(argc == 1) print_crms_usage(argv[0]), exit(EXIT_FAILURE);
    // Build settings struct

    marksplit_settings_t settings = {0};

    settings.hp_threshold = 10;
    settings.n_nucs = DEFAULT_N_NUCS;
    settings.notification_interval = 1000000;
    settings.threads = 1;
    settings.max_blen = -1;
    settings.gzip_compression = 1;
    settings.cleanup = 1;
    sprintf(settings.mode, "wT");

    //omp_set_dynamic(0); // Tell omp that I want to set my number of threads 4realz
    int c;
    while ((c = getopt(argc, argv, "T:t:o:n:s:l:m:r:p:f:v:u:g:i:zwcdh?S=")) > -1) {
        switch(c) {
            case 'c': LOG_WARNING("Deprecated option -c.\n"); break;
            case 'd': settings.run_hash_dmp = 1; break;
            case 'f': settings.ffq_prefix = strdup(optarg); break;
            case 'g': settings.gzip_compression = atoi(optarg); break;
            case 'l': settings.blen = atoi(optarg); break;
            case 'm': settings.offset = atoi(optarg); break;
            case 'n': settings.n_nucs = atoi(optarg); break;
            case 'o': settings.tmp_basename = strdup(optarg); break;
            case 'p': settings.threads = atoi(optarg); break;
            case 'r': settings.rescaler_path = strdup(optarg); break;
            case 's':
                settings.homing_sequence = strdup(optarg);
                settings.homing_sequence_length = strlen(settings.homing_sequence);
                break;
            case 't': settings.hp_threshold = atoi(optarg); break;
            case 'u': settings.notification_interval = atoi(optarg); break;
            case 'v': settings.max_blen = atoi(optarg); break;
            case 'w': settings.cleanup = 0; break;
            case 'z': settings.gzip_output = 1; break;
            case 'T': sprintf(settings.mode, "wb%i", atoi(optarg) % 10); break;
            case 'S': settings.is_se = 1; break;
            case '=': settings.to_stdout = 1; break;
            case '?': case 'h': print_crms_usage(argv[0]), exit(EXIT_SUCCESS);
        }
    }

    if(!settings.gzip_output) settings.gzip_compression = 0;

    // Check for proper command-line usage.
    if(settings.is_se) {
        if(argc < 4) {
            print_crms_usage(argv[0]);
            exit(EXIT_FAILURE);
        }
        if(argc != optind + 1) {
            print_crms_usage(argv[0]);
            LOG_EXIT("Exactly one read fastq required for single-end. See usage.\n");
        }
        // Number of file handles
        settings.n_handles = ipow(4, settings.n_nucs);
        if(settings.n_handles * 2 > get_fileno_limit()) {
            LOG_INFO("Increasing nofile limit from %i to %i.\n", get_fileno_limit(), settings.n_handles * 2);
            increase_nofile_limit(settings.n_handles * 2);
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
            LOG_INFO("Increasing nofile limit from %i to %i.\n", get_fileno_limit(), settings.n_handles * 4);
            increase_nofile_limit(settings.n_handles * 4);
        }
        // Handle filenames
        settings.input_r1_path = strdup(argv[optind]);
        settings.input_r2_path = strdup(argv[optind + 1]);
    }
    // Required parameters
    if(settings.ffq_prefix && !settings.run_hash_dmp)
        LOG_EXIT("Final fastq prefix option provided but run_hash_dmp not selected."
                "Either eliminate the -f flag or add the -d flag.\n");

    // Handle number of threads
    omp_set_num_threads(settings.threads);

    // Handle homing sequence
    if(!settings.homing_sequence)
        LOG_EXIT("Homing sequence not provided. Required.\n");
    clean_homing_sequence(settings.homing_sequence);

    // Handle barcode length
    if(!settings.blen)
        LOG_EXIT("Barcode length not provided. Required. Abort!\n");
    // Handle the offset parameter. If false, blen doesn't change.
    if(settings.is_se) {
        settings.blen -= settings.offset;
        if(settings.max_blen > 0 && settings.max_blen < settings.blen)
            LOG_EXIT("Max blen (%i) must be less than the minimum blen provided (%i).\n",
                    settings.max_blen, settings.blen / 2);
        if(settings.max_blen < 0) settings.max_blen = settings.blen;
    } else {
        settings.blen = (settings.blen - settings.offset) * 2;
        if(settings.max_blen > 0 && settings.max_blen * 2 < settings.blen)
            LOG_EXIT("Max blen (%i) must be less than the minimum blen provided (%i).\n",
                    settings.max_blen, settings.blen / 2);
        settings.blen1_2 = settings.blen / 2;
        if(settings.max_blen < 0) settings.max_blen = settings.blen1_2;
    }

    if(!settings.tmp_basename) {
        // If tmp_basename unset, create a random temporary file prefix.
        kstring_t rs = salted_rand_string(settings.input_r1_path, RANDSTR_SIZE);
        settings.tmp_basename = ks_release(&rs);
        LOG_INFO("Temporary basename not provided. Defaulting to random: %s.\n",
                  settings.tmp_basename);
    }

    // Run core
    mark_splitter_t *splitter = pp_split_inline(&settings);
    //mark_splitter_t *splitter = settings.is_se ? pp_split_inline_se(&settings): pp_split_inline(&settings);
    splitterhash_params_t *params;
    char ffq_r1[500];
    char ffq_r2[500];
    if(!settings.run_hash_dmp) {
        LOG_INFO("mark/split complete.\n");
        goto cleanup;
    }
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
    else cat_fastqs(&settings, params, ffq_r1, ffq_r2);
    cleanup_hashdmp(&settings, params);
    splitterhash_destroy(params);
    cleanup:
    free_marksplit_settings(settings);
    splitter_destroy(splitter);
    LOG_INFO("Successfully completed bmftools dmp!\n");
    return EXIT_SUCCESS;
}

