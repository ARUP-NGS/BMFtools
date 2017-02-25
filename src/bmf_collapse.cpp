#include "bmf_collapse.h"

#include <getopt.h>
#include <omp.h>
#include <zlib.h>
#include "dlib/nix_util.h"
#include "lib/binner.h"
#include "lib/mseq.h"
#define __STDC_FORMAT_MACROS
#include <cinttypes>

namespace bmf {


void dmp_usage() {
    fprintf(stderr, "Collapses initial fastq by exact barcode matching.\n"
                    "Subcommands:\ninline: Inline barcoded chemistry.\n"
                    "secondary: Secondary barcoded chemistry.\n");
}


void idmp_usage()
{
        fprintf(stderr,
                        "Collapses inline barcoded fastq data.\n"
                        "Usage: bmftools collapse inline <options> <r1.fq> <r2.fq>"
                        "\nFlags:\n"
                        "-S: Run in single-end mode. (Ignores read 2)\n"
                        "-=: Emit interleaved final output to stdout.\n"
                        "-l: Number of nucleotides at the beginning of each read to "
                        "use for barcode. Final barcode length is twice this. REQUIRED.\n"
                        "-s: homing sequence. REQUIRED.\n"
                        "-o: Mark/split temporary file basename. Defaults to a random string.\n"
                        "-t: Homopolymer failure threshold. A molecular barcode with"
                        " a homopolymer of length >= this limit is flagged as QC fail."
                        "Default: 10.\n"
                        "-I: Ignore homing sequence. Not recommended, but possible under certain experimental conditions.\n"
                        "-n: Number of nucleotides at the beginning of the barcode to use to split the output. Default: %i.\n"
                        "-m: Mask first n nucleotides in read for barcode. Default: 0.\n"
                        "-p: Number of threads to use if running uthash_dmp. Default: %i.\n"
                        "-D: Use this flag to only mark/split and avoid final demultiplexing/consolidation.\n"
                        "-f: If running hash_dmp, this sets the Final Fastq Prefix. \n"
                        "The Final Fastq files will be named '<ffq_prefix>.R1.fq' and '<ffq_prefix>.R2.fq'.\n"
                        "-r: Path to flat text file with rescaled quality scores. If not provided, it will not be used.\n"
                        "-v: Maximum barcode length for a variable length barcode dataset. If left as default value,"
                        " (-1), other barcode lengths will not be considered.\n"
                        "-z: Flag to write out final output as compressed. Default: False.\n"
                        "-T: If unset, write uncompressed plain text temporary files. If not, use that compression level for temporary files.\n"
                        "-g: Gzip compression ratio if writing gzipped. Default (if writing compressed): 1 (mostly to reduce I/O).\n"
                        "-u: Set notification/update interval for split. Default: 1000000.\n"
                        "-w: Set flag to leave temporary files. Primarily for debugging.\n"
                        "-h: Print usage.\n"
                    , DEFAULT_N_NUCS, DEFAULT_N_THREADS);

}

kstring_t salted_rand_string(char *infname, size_t n_rand) {
    if(strchr(infname, '/')) infname = strrchr(infname, '/') + 1;
    std::string tmp(infname);
    while(strchr(tmp.c_str(), '.')) {
        int n(tmp.c_str() + tmp.size() - strchr(tmp.c_str(), '.') + 1);
        while(n--) tmp.pop_back();
    }
    size_t flen(tmp.size() + n_rand);
    const char cstr[] {"ABCDEFGHIJKLMNOPQRTSUVWXYZ1234567890"};
    while(tmp.size() < flen) tmp.push_back(cstr[rand() % (sizeof(cstr) - 1)]);
    kstring_t ret{0};
    kputs(tmp.data(), &ret);
    return ret;
}
/*
 * Makes an output filename when not provided. Handles cases with and without period in the input fq name.
 */
char *make_salted_fname(char *base)
{
    if(strchr(base, '\0')) {
        kstring_t rs(salted_rand_string(base, RANDSTR_SIZE));
        LOG_INFO("No output final prefix set. Defaulting to variation on input ('%s').\n", base);
        return rs.s;
    } else {
        char *tmp((char *)malloc(RANDSTR_SIZE + 1));
        return dlib::rand_string(tmp, RANDSTR_SIZE);
    }
}

/*
 * Makes an output filename when not provided. Handles cases with and without period in the input fq name.
 */
void make_outfname(marksplit_settings_t *settings)
{
    settings->ffq_prefix = make_salted_fname(settings->input_r1_path);
}

void cleanup_hashdmp(marksplit_settings_t *settings, splitterhash_params_t *params)
{
    if(!settings->cleanup) return;
    #pragma omp parallel for
    for(int i = 0; i < params->n; ++i) {
        kstring_t ks{0, 0, nullptr};
        ksprintf(&ks, "rm %s %s", params->outfnames_r1[i], settings->is_se ? "": params->outfnames_r2[i]);
        dlib::check_call(ks.s);
        free(ks.s);
    }
}

/*
 * Executes hash_dmp_fn on each of the temporary files in the splitterhash
 * and cleans up if not disabled.
 */
void parallel_hash_dmp_core(marksplit_settings_t *settings, splitterhash_params_t *params, hash_dmp_fn func)
{
    #pragma omp parallel for schedule(dynamic, 1)
    for(int i = 0; i < settings->n_handles; ++i) {
        LOG_DEBUG("Now running hash dmp core on input filename %s and output filename %s.\n",
                 params->infnames_r1[i], params->outfnames_r1[i]);
        func(params->infnames_r1[i], params->outfnames_r1[i], settings->gzip_compression);
        if(settings->cleanup) {
            kstring_t ks{0, 0, nullptr};
            ksprintf(&ks, "rm %s", params->infnames_r1[i]);
            dlib::check_call(ks.s);
            free(ks.s);
        }
    }
    if(settings->is_se) {
        LOG_DEBUG("Only dmp'ing single end.\n");
        return; // Don't dmp imaginary files....
    }
    #pragma omp parallel for schedule(dynamic, 1)
    for(int i = 0; i < settings->n_handles; ++i) {
        LOG_DEBUG("Now running hash dmp core on input filename %s and output filename %s.\n",
                  params->infnames_r2[i], params->outfnames_r2[i]);
        func(params->infnames_r2[i], params->outfnames_r2[i], settings->gzip_compression);
        if(settings->cleanup) {
            kstring_t ks{0, 0, nullptr};
            ksprintf(&ks, "rm %s", params->infnames_r2[i]);
            dlib::check_call(ks.s);
            free(ks.s);
        }
    }
}


void cat_fastqs_se(marksplit_settings_t *settings, splitterhash_params_t *params, char *ffq_r1)
{
    kstring_t ks{0, 0, nullptr};
    // Clear output files.
    ksprintf(&ks, settings->gzip_output ? "> %s.gz" : "> %s", ffq_r1);
    dlib::check_call(ks.s);
    ks.l = 0;
    ksprintf(&ks, "/bin/cat ");
    for(int i(0); i < settings->n_handles; ++i) {
        if(!dlib::isfile(params->outfnames_r1[i]))
            LOG_EXIT("Output filename is not a file. Abort! ('%s').\n", params->outfnames_r1[i]);
        ksprintf(&ks, " %s", params->outfnames_r1[i]);
    }
    ksprintf(&ks, " > %s", ffq_r1);
    if(settings->gzip_output) kputsnl(".gz", &ks);
    dlib::check_popen(ks.s);
    free(ks.s);
}
/*
 * Make sure that no rescaler values are invalid
 */
void check_rescaler(marksplit_settings_t *settings, int arr_size)
{
    if(settings->rescaler)
        for(int i(0); i < arr_size; ++i)
            if(settings->rescaler[i] <= 0)
                LOG_EXIT("Invalid value in rescaler %i at index %i.\n", settings->rescaler[i], i);
}
/*
 * Emits final results to stdout
 */

void call_stdout(marksplit_settings_t *settings, splitterhash_params_t *params, char *ffq_r1, char *ffq_r2)
{
    kstring_t str1{0, 0, nullptr}, str2{0, 0, nullptr};
    kstring_t final{0, 0, nullptr};
    kputs((settings->gzip_output)? "zcat": "cat", &str1);
    ks_resize(&str1, 1 << 10);
    for(int i(0); i < settings->n_handles; ++i)
        ksprintf(&str1, " %s", params->outfnames_r1[i]);
    kputsnl(" | paste -d'~' - - - - ", &str1);
    str2.s = dlib::kstrdup(&str1); // strdup the string.
    for(uint32_t i(0); i < str2.l; ++i) {
        LOG_DEBUG("Current str.s + i: %s.\n", str2.s + i);
        if(memcmp(str2.s + i, "R1", 2) == 0)
            str2.s[i + 1] = '2';
    }

    const char final_template[]{"pr -mts'~' <(%s) <(%s) | tr '~' '\\n'"};
    ksprintf(&final, final_template, str1.s, str2.s);
    dlib::bash_system(final.s);
    free(str1.s), free(str2.s);
    free(final.s);
}

void cat_fastqs(marksplit_settings_t *settings, splitterhash_params_t *params, char *ffq_r1, char *ffq_r2)
{
    settings->is_se ? cat_fastqs_se(settings, params, ffq_r1)
                    : cat_fastqs_pe(settings, params, ffq_r1, ffq_r2);
}

void cat_fastqs_pe(marksplit_settings_t *settings, splitterhash_params_t *params, char *ffq_r1, char *ffq_r2)
{
    kstring_t ks1{0, 0, nullptr};
    kputsnl("> ", &ks1), kputs(ffq_r1, &ks1);
    if(settings->gzip_output) kputsnl(".gz", &ks1);
    // Clear output files.
    //ksprintf(&ks1, settings->gzip_output ? "> %s.gz" : "> %s", ffq_r1);
    dlib::check_call(ks1.s); ks1.l = 0;
    ksprintf(&ks1, settings->gzip_output ? "> %s.gz" : "> %s", ffq_r2);
    dlib::check_call(ks1.s); ks1.l = 0;
    kputsnl("/bin/cat ", &ks1);
    kstring_t ks2{0};
    kputsn(ks1.s, ks1.l, &ks2);
    for(int i(0); i < settings->n_handles; ++i) {
        if(!dlib::isfile(params->outfnames_r1[i]))
            LOG_EXIT("Output filename is not a file. Abort! ('%s').\n", params->outfnames_r1[i]);
        if(!dlib::isfile(params->outfnames_r2[i]))
            LOG_EXIT("Output filename is not a file. Abort! ('%s').\n", params->outfnames_r2[i]);
        ksprintf(&ks1, "%s ", params->outfnames_r1[i]);
        ksprintf(&ks2, "%s ", params->outfnames_r2[i]);
    }
    ksprintf(&ks1, " > %s", ffq_r1);
    ksprintf(&ks2, " > %s", ffq_r2);
    if(settings->gzip_output) {
        kputsnl(".gz", &ks1);
        kputsnl(".gz", &ks2);
    }
    FILE *c1_popen(popen(ks1.s, "w"));
    FILE *c2_popen(popen(ks2.s, "w"));
    int ret;
    if((ret = ((pclose(c2_popen) << 8) | pclose(c1_popen)))) {
        LOG_EXIT("Background cat command(s) failed. (Code: %i, '%s' or %i, '%s').\n", ret >> 8, ks1.s, ret & 0xff, ks2.s);
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
        case 'a': case 'c': case 'g': case 't':
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
mark_splitter_t pp_split_inline_se(marksplit_settings_t *settings)
{
    uint64_t bin(0), count(1); // Initialze count to 1 because of the read used to get read length.
    const int default_nlen(settings->blen + settings->offset + settings->homing_sequence_length);
    LOG_DEBUG("Opening fastq file %s.\n", settings->input_r1_path);
    if(!dlib::isfile(settings->input_r1_path))
        LOG_EXIT("Could not open read paths: at least one is not a file.\n");
    if(settings->rescaler_path)
        settings->rescaler = parse_1d_rescaler(settings->rescaler_path);
    mark_splitter_t splitter(init_splitter(settings));
    gzFile fp(gzopen(settings->input_r1_path, "r"));
    kseq_t *seq(kseq_init(fp));
    // Manually process the first pair of reads so that we have the read length.
    int pass_fail = 1;
    if(kseq_read(seq) < 0) {
            free_marksplit_settings(*settings);
            splitter_destroy(&splitter);
            LOG_EXIT("Could not open fastqs for reading. Abort!\n");
    }
    LOG_DEBUG("Read length (inferred): %lu.\n", seq->seq.l);
    check_rescaler(settings, seq->seq.l * 4 * 2 * NQSCORES);
    tmp_mseq_t *tmp(init_tm_ptr(seq->seq.l, settings->blen));
    int n_len(nlen_homing_se(seq, settings, default_nlen, &pass_fail));
    mseq_t *rseq(mseq_rescale_init(seq, settings->rescaler, tmp, 0));
    rseq->barcode[settings->blen] = '\0';
    memcpy(rseq->barcode, seq->seq.s + settings->offset, settings->blen);
    pass_fail &= test_hp(rseq->barcode, settings->hp_threshold);
    // Get first barcode.
    update_mseq(rseq, seq, settings->rescaler, tmp, n_len, 0);
    bin = get_binner_type(rseq->barcode, settings->n_nucs, uint64_t);
    assert(bin < (uint64_t)settings->n_handles);
    mseq2fq_stranded(splitter.tmp_out_handles_r1[bin], rseq, pass_fail, rseq->barcode, 'F');
    while(LIKELY(kseq_read(seq) >= 0)) {
        if(UNLIKELY(++count % settings->notification_interval == 0))
            LOG_INFO("Number of records processed: %lu.\n", count);
        // Sets pass_fail and gets n_len
        n_len = nlen_homing_se(seq, settings, default_nlen, &pass_fail);
        // Update mseq
        update_mseq(rseq, seq, settings->rescaler, tmp, n_len, 0);
        // Update barcode
        memcpy(rseq->barcode, seq->seq.s + settings->offset, settings->blen);
        // Update QC Fail
        pass_fail &= test_hp(rseq->barcode, settings->hp_threshold);
        // Get bin
        bin = bmf::get_binner_type(rseq->barcode, settings->n_nucs, uint64_t);
        assert(bin < (uint64_t)settings->n_handles);
        // Write the processed read to the bin
        mseq2fq_stranded(splitter.tmp_out_handles_r1[bin], rseq, pass_fail, rseq->barcode, 'F');
    }
    LOG_INFO("Collapsing %lu initial reads....\n", count);
    LOG_DEBUG("Cleaning up.\n");
    for(int i(0); i < splitter.n_handles; ++i) gzclose(splitter.tmp_out_handles_r1[i]);
    tm_destroy(tmp);
    mseq_destroy(rseq);
    kseq_destroy(seq);
    gzclose(fp);
    return splitter;
}


/*
 * Pre-processes (pp) and splits fastqs with inline barcodes.
 */
mark_splitter_t pp_split_inline(marksplit_settings_t *settings)
{
    LOG_DEBUG("Opening fastq files %s and %s.\n", settings->input_r1_path, settings->input_r2_path);
    if(!(strcmp(settings->input_r1_path, settings->input_r2_path))) {
        LOG_EXIT("Hey, it looks like you're trying to use the same path for both r1 and r2. "
                "At least try to fool me by making a symbolic link.\n");
    }
    if(!dlib::isfile(settings->input_r1_path) || !dlib::isfile(settings->input_r2_path)) {
        LOG_EXIT("Could not open read paths: at least one is not a file.\n");
    }
    if(settings->rescaler_path) settings->rescaler = parse_1d_rescaler(settings->rescaler_path);
    mark_splitter_t splitter(init_splitter(settings));
    gzFile fp1(gzopen(settings->input_r1_path, "r"));
    gzFile fp2(gzopen(settings->input_r2_path, "r"));
    kseq_t *seq1(kseq_init(fp1));
    kseq_t *seq2(kseq_init(fp2));
    // Manually process the first pair of reads so that we have the read length.
    int pass_fail(1);
    if(kseq_read(seq1) < 0 || kseq_read(seq2) < 0) {
            free_marksplit_settings(*settings);
            splitter_destroy(&splitter);
            LOG_EXIT("Could not open fastqs for reading. Abort!\n");
    }
    LOG_DEBUG("Read length (inferred): %lu.\n", seq1->seq.l);
    check_rescaler(settings, seq1->seq.l * 4 * 2 * NQSCORES);
    tmp_mseq_t *tmp(init_tm_ptr(seq1->seq.l, settings->blen));
    int switch_reads(switch_test(seq1, seq2, settings->offset));
    const int default_nlen(settings->blen1_2 + settings->offset + settings->homing_sequence_length);
    int n_len(settings->ignore_homing ? settings->blen1_2 + settings->offset + settings->homing_sequence_length
                                      : nlen_homing_default(seq1, seq2, settings, default_nlen, &pass_fail));
    mseq_t *rseq1(mseq_rescale_init(seq1, settings->rescaler, tmp, 0));
    mseq_t *rseq2(mseq_rescale_init(seq2, settings->rescaler, tmp, 1));
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
    uint64_t bin(get_binner_type(rseq1->barcode, settings->n_nucs, uint64_t));
    assert(bin < (uint64_t)settings->n_handles);
    if(switch_reads) {
        mseq2fq_stranded(splitter.tmp_out_handles_r1[bin], rseq2, pass_fail, rseq1->barcode, 'R');
        mseq2fq_stranded(splitter.tmp_out_handles_r2[bin], rseq1, pass_fail, rseq1->barcode, 'R');
    } else {
        mseq2fq_stranded(splitter.tmp_out_handles_r1[bin], rseq1, pass_fail, rseq1->barcode, 'F');
        mseq2fq_stranded(splitter.tmp_out_handles_r2[bin], rseq2, pass_fail, rseq1->barcode, 'F');
    }
    uint64_t count(1uL);
    while(LIKELY(kseq_read(seq1) >= 0 && kseq_read(seq2) >= 0)) {
        if(UNLIKELY(++count % settings->notification_interval == 0))
            LOG_INFO("Number of records processed: %lu.\n", count);
        // Sets pass_fail
        n_len = settings->ignore_homing ? settings->blen1_2 + settings->offset
                                        : nlen_homing_default(seq1, seq2, settings, default_nlen, &pass_fail);
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
            mseq2fq_stranded(splitter.tmp_out_handles_r1[bin], rseq2, pass_fail, rseq1->barcode, 'R');
            mseq2fq_stranded(splitter.tmp_out_handles_r2[bin], rseq1, pass_fail, rseq1->barcode, 'R');
        } else {
            memcpy(rseq1->barcode, seq1->seq.s + settings->offset, settings->blen1_2);
            memcpy(rseq1->barcode + settings->blen1_2, seq2->seq.s + settings->offset, settings->blen1_2);
            pass_fail &= test_hp(rseq1->barcode, settings->hp_threshold);
            bin = bmf::get_binner_type(rseq1->barcode, settings->n_nucs, uint64_t);
            assert(bin < (uint64_t)settings->n_handles);
            mseq2fq_stranded(splitter.tmp_out_handles_r1[bin], rseq1, pass_fail, rseq1->barcode, 'F');
            mseq2fq_stranded(splitter.tmp_out_handles_r2[bin], rseq2, pass_fail, rseq1->barcode, 'F');
        }
    }
    LOG_INFO("Collapsing %lu initial read pairs....\n", count);
    LOG_DEBUG("Cleaning up.\n");
    for(int i(0); i < splitter.n_handles; ++i) {
        gzclose(splitter.tmp_out_handles_r1[i]);
        gzclose(splitter.tmp_out_handles_r2[i]);
    }
    tm_destroy(tmp);
    mseq_destroy(rseq1), mseq_destroy(rseq2);
    kseq_destroy(seq1), kseq_destroy(seq2);
    gzclose(fp1), gzclose(fp2);
    return splitter;
}

int idmp_main(int argc, char *argv[]);
extern int sdmp_main(int argc, char *argv[]);

int collapse_main(int argc, char *argv[])
{
    if(argc == 1) dmp_usage(), exit(EXIT_FAILURE);
    if(strcmp(argv[1], "inline") == 0)
        return idmp_main(argc - 1, argv + 1);
    if(strcmp(argv[1], "secondary") == 0)
        return sdmp_main(argc - 1, argv + 1);
    fprintf(stderr, "Unrecognized bmftools dmp subcommand %s. Abort!\n", argv[1]);
    return EXIT_FAILURE;
}


int idmp_main(int argc, char *argv[])
{
    if(argc == 1) idmp_usage(), exit(EXIT_FAILURE);
    // Build settings struct

    marksplit_settings_t settings{0};

    settings.hp_threshold = 10;
    settings.n_nucs = DEFAULT_N_NUCS;
    settings.notification_interval = 1000000;
    settings.threads = DEFAULT_N_THREADS;
    settings.max_blen = -1;
    settings.gzip_compression = 1;
    settings.cleanup = 1;
    settings.run_hash_dmp = 1;
#if ZLIB_VER_MAJOR <= 1 && ZLIB_VER_MINOR <= 2 && ZLIB_VER_REVISION < 5
#pragma message("Note: zlib version < 1.2.5 doesn't support transparent file writing. Writing uncompressed temporary gzip files by default.")
    // If not set, zlib compresses all our files enormously.
    sprintf(settings.mode, "wb0");
#else
    sprintf(settings.mode, "wT");
#endif

    //omp_set_dynamic(0); // Tell omp that I want to set my number of threads 4realz
    int c;
    while ((c = getopt(argc, argv, "T:t:o:n:s:l:m:r:p:f:v:u:g:i:zwcdDh?S=")) > -1) {
        switch(c) {
            case 'c': LOG_WARNING("Deprecated option -c.\n"); break;
            case 'd': LOG_WARNING("Deprecated option -d.\n"); break;
            case 'D': settings.run_hash_dmp = 0; break;
            case 'f': settings.ffq_prefix = strdup(optarg); break;
            case 'g': settings.gzip_compression = (uint32_t)atoi(optarg)%10; break;
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
            case '?': case 'h': idmp_usage(); exit(EXIT_SUCCESS);
        }
    }

    if(!settings.gzip_output) settings.gzip_compression = 0;

    // Check for proper command-line usage.
    if(settings.is_se) {
        if(argc < 4) {
            idmp_usage();
            exit(EXIT_FAILURE);
        }
        if(argc < optind + 1) {
            idmp_usage();
        }
        if(argc == optind + 2) {
            LOG_WARNING("Note: two read paths were provided but single-end mode was selected.\n");
        }
        // Number of file handles
        settings.n_handles = dlib::ipow(4, settings.n_nucs);
        if(settings.n_handles * 2 > dlib::get_fileno_limit()) {
            LOG_INFO("Increasing nofile limit from %i to %i.\n", dlib::get_fileno_limit(), settings.n_handles * 2);
            dlib::increase_nofile_limit(settings.n_handles * 2);
        }
        // Handle filenames
        settings.input_r1_path = strdup(argv[optind]);
    } else {
        if(argc < 5) idmp_usage(), exit(EXIT_FAILURE);
        if(argc != optind + 2) {
            fprintf(stderr, "[E:%s] Both read 1 and read 2 fastqs are required for paired-end. See usage.\n", __func__);
            idmp_usage();
            return EXIT_FAILURE;
        }
        // Number of file handles
        settings.n_handles = dlib::ipow(4, settings.n_nucs);
        if(settings.n_handles * 4 > dlib::get_fileno_limit()) {
            LOG_INFO("Increasing nofile limit from %i to %i.\n", dlib::get_fileno_limit(), settings.n_handles * 4);
            dlib::increase_nofile_limit(settings.n_handles * 4);
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
    if(!settings.homing_sequence && !settings.ignore_homing)
        LOG_EXIT("Homing sequence not provided. Required.\n");
    if(settings.homing_sequence) clean_homing_sequence(settings.homing_sequence);

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
        kstring_t rs(salted_rand_string(settings.input_r1_path, RANDSTR_SIZE));
        settings.tmp_basename = ks_release(&rs);
        LOG_INFO("Temporary basename not provided. Defaulting to random: %s.\n",
                  settings.tmp_basename);
    }

    // Run core
    mark_splitter_t splitter(settings.is_se ? pp_split_inline_se(&settings)
                                            : pp_split_inline(&settings));
    splitterhash_params_t *params(init_splitterhash(&settings, &splitter));
    kstring_t ffq_r1{0,0,nullptr};
    kstring_t ffq_r2{0,0,nullptr};
    if(!settings.run_hash_dmp) {
        LOG_INFO("mark/split complete.\n");
        goto cleanup;
    }
    if(!settings.ffq_prefix) make_outfname(&settings);
    // Run cores.
    parallel_hash_dmp_core(&settings, params, &stranded_hash_dmp_core);

    // Remove temporary split files.
    ksprintf(&ffq_r1, "%s.R1.fq", settings.ffq_prefix);
    ksprintf(&ffq_r2, "%s.R2.fq", settings.ffq_prefix);
    // Cat temporary files together.
    if(settings.to_stdout)
        call_stdout(&settings, params, ffq_r1.s, ffq_r2.s);
    else cat_fastqs(&settings, params, ffq_r1.s, ffq_r2.s);
    free(ffq_r1.s), free(ffq_r2.s);
    cleanup_hashdmp(&settings, params);
    splitterhash_destroy(params);
    cleanup:
    free_marksplit_settings(settings);
    splitter_destroy(&splitter);
    LOG_INFO("Successfully completed bmftools collapse inline!\n");
    return EXIT_SUCCESS;
} /* idmp_main */


void sdmp_usage(char *argv[])
{
        fprintf(stderr,
                        "Performs molecular demultiplexing for secondary index barcoded experiments.\n"
                        "Usage: bmftools collapse secondary <options> <r1.fq> <r2.fq>\n"
                        "Flags:\n"
                        "-i: Index fastq path. REQUIRED.\n"
                        "-t: Homopolymer failure threshold. A molecular barcode with a homopolymer of length >= this limit is flagged as QC fail. Default: 10\n"
                        "-o: Temporary fastq file prefix.\n"
                        "-n: Number of nucleotides at the beginning of the barcode to use to split the output. Default: %i.\n"
                        "-z: Flag to write gzip compressed output. Default: False.\n"
                        "-T: If unset, write uncompressed plain text temporary files. If not, use that compression level for temporary files.\n"
                        "-g: Gzip compression ratio if writing compressed. Default: 1 (mostly to reduce I/O).\n"
                        "-s: Number of bases from reads 1 and 2 with which to salt the barcode. Default: 0.\n"
                        "-m: Number of bases in the start of reads to skip when salting. Default: 0.\n"
                        "-D: Use this flag to only mark/split and avoid final demultiplexing/consolidation.\n"
                        "-p: Number of threads to use if running hash_dmp. Default: %i.\n"
                        "-v: Set notification interval for split. Default: 1000000.\n"
                        "-r: Path to flat text file with rescaled quality scores. If not provided, it will not be used.\n"
                        "-w: Flag to leave temporary files instead of deleting them, as in default behavior.\n"
                        "-f: If running hash_dmp, this sets the Final Fastq Prefix. \n"
                        "-S: Single-end mode. Ignores read 2.\n"
                        "-=: Emit final fastqs to stdout in interleaved form. Ignores -f.\n"
                , DEFAULT_N_NUCS, DEFAULT_N_THREADS);
}

static mark_splitter_t splitmark_core_rescale(marksplit_settings_t *settings)
{
    LOG_DEBUG("Path to index fq: %s.\n", settings->index_fq_path);
    gzFile fp_read1, fp_read2, fp_index;
    kseq_t *seq1, *seq2, *seq_index;
    int l1, l2, l_index;
    mark_splitter_t splitter(init_splitter(settings));
    for(const auto path: {settings->input_r1_path, settings->input_r2_path, settings->index_fq_path})
        if(!dlib::isfile(path))
            LOG_EXIT("%s is not a file. Abort!\n", path);
    // Open fastqs
    LOG_DEBUG("Splitter now opening files R1 ('%s'), R2 ('%s'), index ('%s').\n",
              settings->input_r1_path, settings->input_r2_path, settings->index_fq_path);
    fp_read1 = gzopen(settings->input_r1_path, "r"), fp_read2 = gzopen(settings->input_r2_path, "r");
    seq1 = kseq_init(fp_read1), seq2 = kseq_init(fp_read2);
    l1 = kseq_read(seq1), l2 = kseq_read(seq2);

    fp_index = gzopen(settings->index_fq_path, "r");
    seq_index = kseq_init(fp_index),
    l_index = kseq_read(seq_index);

    uint64_t bin(0);
    int pass_fail(1);
    tmp_mseq_t *tmp(init_tm_ptr(seq1->seq.l, seq_index->seq.l + 2 * settings->salt));
    if(l1 < 0 || l2 < 0 || l_index < 0)
        LOG_EXIT("Could not read input fastqs. Abort mission!\n");
    check_rescaler(settings, NQSCORES * 2 * 4 * seq1->seq.l);
    mseq_t *rseq1(mseq_init(seq1, settings->rescaler, 0)); // rseq1 is initialized
    mseq_t *rseq2(mseq_init(seq2, settings->rescaler, 1)); // rseq2 is initialized
    memcpy(rseq1->barcode, seq1->seq.s + settings->offset, settings->salt); // Copy in the appropriate nucleotides.
    memcpy(rseq1->barcode + settings->salt, seq_index->seq.s, seq_index->seq.l); // Copy in the barcode
    memcpy(rseq1->barcode + settings->salt + seq_index->seq.l, seq2->seq.s + settings->offset, settings->salt);
    rseq1->barcode[settings->salt * 2 + seq_index->seq.l] = '\0';
    update_mseq(rseq1, seq1, settings->rescaler, tmp, 0, 0);
    update_mseq(rseq2, seq2, settings->rescaler, tmp, 0, 1);
    pass_fail = test_hp(rseq1->barcode, settings->hp_threshold);
    bin = get_binner_type(rseq1->barcode, settings->n_nucs, uint64_t);
    mseq2fq(splitter.tmp_out_handles_r1[bin], rseq1, pass_fail, rseq1->barcode);
    mseq2fq(splitter.tmp_out_handles_r2[bin], rseq2, pass_fail, rseq1->barcode);
    uint64_t count(1uL);
    while(LIKELY((l1 = kseq_read(seq1)) >= 0 && (l2 = kseq_read(seq2) >= 0))
            && (l_index = kseq_read(seq_index)) >= 0) {
        if(UNLIKELY(++count % settings->notification_interval == 0))
            LOG_INFO("Number of records processed: %lu.\n", count);
        memcpy(rseq1->barcode, seq1->seq.s + settings->offset, settings->salt); // Copy in the appropriate nucleotides.
        memcpy(rseq1->barcode + settings->salt, seq_index->seq.s, seq_index->seq.l); // Copy in the barcode
        memcpy(rseq1->barcode + settings->salt + seq_index->seq.l, seq2->seq.s + settings->offset, settings->salt);
        update_mseq(rseq1, seq1, settings->rescaler, tmp, 0, 0);
        update_mseq(rseq2, seq2, settings->rescaler, tmp, 0, 1);
        pass_fail = test_hp(rseq1->barcode, settings->hp_threshold);
        bin = get_binner_type(rseq1->barcode, settings->n_nucs, uint64_t);
        mseq2fq(splitter.tmp_out_handles_r1[bin], rseq1, pass_fail, rseq1->barcode);
        mseq2fq(splitter.tmp_out_handles_r2[bin], rseq2, pass_fail, rseq1->barcode);
    }
    tm_destroy(tmp);
    mseq_destroy(rseq1); mseq_destroy(rseq2);
    kseq_destroy(seq1); kseq_destroy(seq2); kseq_destroy(seq_index);
    gzclose(fp_read1); gzclose(fp_read2); gzclose(fp_index);
    for(int j(0); j < settings->n_handles; ++j) {
        gzclose(splitter.tmp_out_handles_r1[j]);
        gzclose(splitter.tmp_out_handles_r2[j]);
        splitter.tmp_out_handles_r1[j] = splitter.tmp_out_handles_r2[j] = nullptr;
    }
    LOG_INFO("Collapsing %lu initial read pairs....\n", count);
    return splitter;
}

static mark_splitter_t splitmark_core_rescale_se(marksplit_settings_t *settings)
{
    mark_splitter_t splitter(init_splitter(settings));
    if(!dlib::isfile(settings->input_r1_path) ||
       !dlib::isfile(settings->index_fq_path)) {
        LOG_EXIT("At least one input path ('%s', '%s') is not a file. Abort!\n",
                 settings->input_r1_path, settings->index_fq_path);
    }
    // Open fastqs
    gzFile fp(gzopen(settings->input_r1_path, "r")), fp_index(gzopen(settings->index_fq_path, "r"));
    kseq_t *seq(kseq_init(fp)), *seq_index(kseq_init(fp_index));
    int l(kseq_read(seq));
    int l_index(kseq_read(seq_index));

    tmp_mseq_t *tmp(init_tm_ptr(seq->seq.l, seq_index->seq.l + settings->salt));
    if(l < 0 || l_index < 0)
        LOG_EXIT("Could not read input fastqs. Abort mission!\n");
    mseq_t *rseq(mseq_init(seq, settings->rescaler, 0));
    memcpy(rseq->barcode, seq->seq.s + settings->offset, settings->salt); // Copy in the appropriate nucleotides.
    memcpy(rseq->barcode + settings->salt, seq_index->seq.s, seq_index->seq.l); // Copy in the barcode
    rseq->barcode[settings->salt + seq_index->seq.l] = '\0';
    update_mseq(rseq, seq, settings->rescaler, tmp, 0, 0);
    mseq2fq(splitter.tmp_out_handles_r1[get_binner_type(rseq->barcode, settings->n_nucs, uint64_t)],
            rseq, test_hp(rseq->barcode, settings->hp_threshold), rseq->barcode);
    uint64_t count(1uL);
    while (LIKELY((l = kseq_read(seq)) >= 0 && (l_index = kseq_read(seq_index)) >= 0)) {
        if(UNLIKELY(++count % settings->notification_interval == 0))
            fprintf(stderr, "[%s] Number of records processed: %" PRIu64 ".\n", __func__, count);
        memcpy(rseq->barcode, seq->seq.s + settings->offset, settings->salt); // Copy in the appropriate nucleotides.
        memcpy(rseq->barcode + settings->salt, seq_index->seq.s, seq_index->seq.l); // Copy in the barcode
        update_mseq(rseq, seq, settings->rescaler, tmp, 0, 0);
        mseq2fq(splitter.tmp_out_handles_r1[get_binner_type(rseq->barcode, settings->n_nucs, uint64_t)],
                rseq, test_hp(rseq->barcode, settings->hp_threshold), rseq->barcode);
    }
    tm_destroy(tmp);
    mseq_destroy(rseq);
    kseq_destroy(seq); kseq_destroy(seq_index);
    gzclose(fp); gzclose(fp_index);
    for(int j(0); j < settings->n_handles; ++j) {
        gzclose(splitter.tmp_out_handles_r1[j]);
        splitter.tmp_out_handles_r1[j] = nullptr;
    }
    LOG_INFO("Collapsing %lu initial reads....\n", count);
    // Set out handles to nullptr.
    return splitter;
}

int sdmp_main(int argc, char *argv[])
{
    if(argc < 3 || strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0) {
        sdmp_usage(argv); return EXIT_FAILURE;
    }
    // Build settings struct
    marksplit_settings_t settings{0};
    settings.hp_threshold = 10;
    settings.n_nucs = DEFAULT_N_NUCS;
    settings.notification_interval = 1000000;
    settings.threads = DEFAULT_N_THREADS;
    settings.gzip_compression = 1;
    settings.cleanup = 1;
    settings.run_hash_dmp = 1;
#if ZLIB_VER_MAJOR <= 1 && ZLIB_VER_MINOR <= 2 && ZLIB_VER_REVISION < 5
#pragma message("Note: zlib version < 1.2.5 doesn't support transparent file writing. Writing uncompressed temporary gzip files by default.")
    // If not set, zlib compresses all our files enormously.
    sprintf(settings.mode, "wb0");
#else
    sprintf(settings.mode, "wT");
#endif

    int c;
    while ((c = getopt(argc, argv, "t:o:i:n:m:s:f:u:p:g:v:r:T:IhdDczw?S=")) > -1) {
        switch(c) {
            case 'd': LOG_WARNING("Deprecated option -d.\n"); break;
            case 'D': settings.run_hash_dmp = 0; break;
            case 'f': settings.ffq_prefix = strdup(optarg); break;
            case 'i': settings.index_fq_path = strdup(optarg); break;
            case 'I': settings.ignore_homing = 1; break;
            case 'm': settings.offset = atoi(optarg); break;
            case 'n': settings.n_nucs = atoi(optarg); break;
            case 'o': settings.tmp_basename = strdup(optarg);break;
            case 'T': sprintf(settings.mode, "wb%i", atoi(optarg) % 10); break;
            case 'p': settings.threads = atoi(optarg); break;
            case 's': settings.salt = atoi(optarg); break;
            case 't': settings.hp_threshold = atoi(optarg); break;
            case 'v': settings.notification_interval = atoi(optarg); break;
            case 'z': settings.gzip_output = 1; break;
            case 'g': settings.gzip_compression = (uint32_t)atoi(optarg)%10; break;
            case 'w': settings.cleanup = 0; break;
            case 'r':
                settings.rescaler_path = strdup(optarg);
                settings.rescaler = parse_1d_rescaler(settings.rescaler_path);
                break;
            case 'S': settings.is_se = 1; break;
            case '=': settings.to_stdout = 1; break;
            case '?': case 'h': sdmp_usage(argv); return EXIT_SUCCESS;
        }
    }

    dlib::increase_nofile_limit(settings.threads);
    omp_set_num_threads(settings.threads);

    settings.n_handles = dlib::ipow(4, settings.n_nucs);
    if(settings.n_handles * 3 > dlib::get_fileno_limit()) {
        const int o_fnl(dlib::get_fileno_limit());
        dlib::increase_nofile_limit(kroundup32(settings.n_handles));
        fprintf(stderr, "Increased nofile limit from %i to %i.\n", o_fnl,
                kroundup32(settings.n_handles));
    }

    if(argc == 1) {
        sdmp_usage(argv);
        return EXIT_SUCCESS;
    }

    if(settings.is_se) {
        if(argc != optind + 1) {
            fprintf(stderr, "[E:%s] Precisely one input fastq required for se mode. See usage.\n", __func__);
            sdmp_usage(argv);
            return EXIT_FAILURE;
        }
        settings.input_r1_path =  strdup(argv[optind]);
    }
    else {
        if(argc != optind + 2) {
            fprintf(stderr, "[E:%s] Both read 1 and read 2 fastqs are required. See usage.\n", __func__);
            sdmp_usage(argv);
            return EXIT_FAILURE;
        }
        settings.input_r1_path =  strdup(argv[optind]);
        settings.input_r2_path =  strdup(argv[optind + 1]);
    }

    if(!settings.index_fq_path) {
        fprintf(stderr, "[E:%s] Index fastq required. See usage.\n", __func__);
        sdmp_usage(argv);
        return EXIT_FAILURE;
    }
    if(!settings.tmp_basename) {
        settings.tmp_basename = make_salted_fname(settings.input_r1_path);
        fprintf(stderr, "[%s] Mark/split prefix not provided. Defaulting to random string ('%s').\n",
                __func__, settings.tmp_basename);
    }

    splitterhash_params_t *params(nullptr);
    mark_splitter_t splitter(settings.is_se ? splitmark_core_rescale_se(&settings)
                                            : splitmark_core_rescale(&settings));
    if(!settings.run_hash_dmp) {
        fprintf(stderr, "[%s] Finished mark/split. Skipping dmp.\n", __func__);
        goto cleanup;
    }
    if(!settings.ffq_prefix) make_outfname(&settings);
    params = init_splitterhash(&settings, &splitter);
    fprintf(stderr, "[%s] Running dmp block in parallel with %i threads.\n", __func__, settings.threads);

    parallel_hash_dmp_core(&settings, params, &hash_dmp_core);
    // Make sure that both files are empty.
    char ffq_r1[200], ffq_r2[200];
    sprintf(ffq_r1, settings.gzip_output ? "%s.R1.fq": "%s.R1.fq.gz", settings.ffq_prefix);
    sprintf(ffq_r2, settings.gzip_output ? "%s.R2.fq": "%s.R2.fq.gz", settings.ffq_prefix);
    if(settings.to_stdout)
        call_stdout(&settings, params, ffq_r1, ffq_r2);
    else cat_fastqs(&settings, params, ffq_r1, ffq_r2);
    cleanup_hashdmp(&settings, params);
    splitterhash_destroy(params);

    cleanup:
    splitter_destroy(&splitter);
    free_marksplit_settings(settings);
    LOG_INFO("Successfully completed bmftools collapse secondary!\n");
    return EXIT_SUCCESS;
} /* sdmp_main */

} /* namespace bmf */
