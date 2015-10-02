/*
 * Conditional reverse complement Rescale Mark Split
 * For inline Loeb adapters only.
 */

#include <zlib.h>
#include <stdio.h>
#include <stdint.h>
#include <regex.h>
#include <getopt.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include <stdlib.h>
#include <sys/resource.h>
#include "dmp_interface.h"
#include "include/array_parser.h"
#include "include/nix_resource.h"
#include "crms.h"
#include "uthash_dmp_core.c"

#define CAT_BUFFER_SIZE 2000000


void print_crms_usage(char *argv[])
{
        fprintf(stderr, "Usage: %s <options> <Fq.R1.seq> <Fq.R2.seq>"
                        "\nFlags:\n"
                        "-l: Number of nucleotides at the beginning of each read to "
                        "use for barcode. Final barcode length is twice this. REQUIRED.\n"
                        "-s: homing sequence. REQUIRED.\n"
                        "-o: Output basename. Defaults to a variation on input filename.\n"
                        "-t: Homopolymer failure threshold. A molecular barcode with"
                        " a homopolymer of length >= this limit is flagged as QC fail."
                        "Default: 10.\n"
                        "-n: Number of nucleotides at the beginning of the barcode to use to split the output. Default: 4.\n"
                        "-m: Mask first n nucleotides in read for barcode. Default: 0. Recommended: 1.\n"
                        "-p: Number of threads to use if running uthash_dmp.\n"
                        "-d: Use this flag to to run hash_dmp.\n"
                        "-f: If running hash_dmp, this sets the Final Fastq Prefix. \n"
                        "The Final Fastq files will be named '<ffq_prefix>.R1.fq' and '<ffq_prefix>.R2.fq'.\n"
                        "-h: Print usage.\n", argv[0], argv[0]);
}

void print_crms_opt_err(char *argv[], char *optarg)
{
    print_crms_usage(argv);
    fprintf(stderr, "Unrecognized option %s. Abort!\n", optarg);
    exit(1);
}


mark_splitter_t init_splitter_inline(mssi_settings_t* settings_ptr)
{
    mark_splitter_t ret = {
        .n_handles = ipow(4, settings_ptr->n_nucs),
        .n_nucs = settings_ptr->n_nucs,
        .fnames_r1 = NULL,
        .fnames_r2 = NULL,
        .tmp_out_handles_r1 = NULL,
        .tmp_out_handles_r2 = NULL
    };
    // Avoid passing more variables.
    ret.tmp_out_handles_r1 = (FILE **)malloc(settings_ptr->n_handles * sizeof(FILE *));
    ret.tmp_out_handles_r2 = (FILE **)malloc(settings_ptr->n_handles * sizeof(FILE *));
    ret.fnames_r1 = (char **)malloc(ret.n_handles * sizeof(char *));
    ret.fnames_r2 = (char **)malloc(ret.n_handles * sizeof(char *));
    char tmp_buffer [METASYNTACTIC_FNAME_BUFLEN];
    for (int i = 0; i < ret.n_handles; i++) {
        sprintf(tmp_buffer, "%s.tmp.%i.R1.fastq", settings_ptr->output_basename, i);
        ret.fnames_r1[i] = strdup(tmp_buffer);
        sprintf(tmp_buffer, "%s.tmp.%i.R2.fastq", settings_ptr->output_basename, i);
        ret.fnames_r2[i] = strdup(tmp_buffer);
        ret.tmp_out_handles_r1[i] = fopen(ret.fnames_r1[i], "w");
        ret.tmp_out_handles_r2[i] = fopen(ret.fnames_r2[i], "w");
    }
    return ret;
}


/*
 * Pre-processes (pp) and splits fastqs with inline barcodes.
void pp_split_inline(kseq_t *seq1, kseq_t *seq2,
                     mssi_settings_t settings, mark_splitter_t splitter)
{
    int l1, l2;
    uint64_t bin;
#if LOGGING
    int count = 0;
#endif
    char pass_fail;
    int readlen = 0;
    char * barcode;
    int blen1_2 = settings.blen / 2;
    int n_len = blen1_2 + settings.homing_sequence_length + settings.offset;
    barcode = (char *)malloc((settings.blen + 1) * sizeof(char));
    barcode[settings.blen] = '\0'; // Null-terminate
    l1 = kseq_read(seq1);
    l2 = kseq_read(seq2);
    if(l1 < 0 || l2 < 0) {
            fprintf(stderr, "Could not open fastqs for reading. Abort!\n");
            FREE_MSSI_SETTINGS(settings);
            FREE_SPLITTER(splitter);
            exit(EXIT_FAILURE);
    }
    readlen = seq1->seq.l;
    tmp_mseq_t tmp = init_tmp_mseq(readlen, settings.blen);
    // Get first barcode.
    set_barcode(seq1, seq2, barcode, settings.offset, blen1_2);
    pass_fail = test_homing_seq(seq1, seq2, &settings) ? test_hp_inline(barcode, settings.blen, settings.hp_threshold) : '0';
    mseq_t mvar1 = init_rescale_revcmp_mseq(seq1, barcode, settings.rescaler, &tmp, n_len, 0);
    mseq_t mvar2 = init_rescale_revcmp_mseq(seq2, barcode, settings.rescaler, &tmp, n_len, 1);
    bin = get_binnerul(barcode, settings.n_nucs);
    mseq2fq_inline(splitter.tmp_out_handles_r1[bin], &mvar1, pass_fail);
    mseq2fq_inline(splitter.tmp_out_handles_r2[bin], &mvar2, pass_fail);
    do {
#if LOGGING
        count += 1;
        if(!(count % settings.notification_interval)) fprintf(stderr, "Number of records processed: %i.\n", count);
#endif
        // Iterate through second fastq file.
        set_barcode(seq1, seq2, barcode, settings.offset, blen1_2);
        pass_fail = test_homing_seq(seq1, seq2, &settings) ? test_hp_inline(barcode, settings.blen, settings.hp_threshold) : '0';
        //fprintf(stdout, "Randomly testing to see if the reading is working. %s", seq1->seq.s);
        update_mseq(&mvar1, barcode, seq1, settings.rescaler, &tmp, n_len, 0);
        update_mseq(&mvar2, barcode, seq2, settings.rescaler, &tmp, n_len, 1);
        bin = get_binnerul(barcode, settings.n_nucs);
        mseq2fq_inline(splitter.tmp_out_handles_r1[bin], &mvar1, pass_fail);
        mseq2fq_inline(splitter.tmp_out_handles_r2[bin], &mvar2, pass_fail);
    } while (((l1 = kseq_read(seq1)) >= 0) && ((l2 = kseq_read(seq2)) >= 0));
    for(int i = 0; i < splitter.n_handles; ++i) {
        fclose(splitter.tmp_out_handles_r1[i]);
        fclose(splitter.tmp_out_handles_r2[i]);
    }
    tmp_mseq_destroy(tmp);
    mseq_destroy(&mvar1);
    mseq_destroy(&mvar2);
    free(barcode);
    return;
}
 */


/*
 * Pre-processes (pp) and splits fastqs with inline barcodes.
 */
mark_splitter_t *pp_split_inline1(char *r1fq, char *r2fq,
                     mssi_settings_t *settings)
{
    fprintf(stderr, "Now beginning pp_split_inline with fastq paths %s and %s.\n", r1fq, r2fq);
    mark_splitter_t *splitter = (mark_splitter_t *)malloc(sizeof(mark_splitter_t));
    *splitter = init_splitter_inline(settings);
    gzFile fp1 = gzopen(r1fq, "r");
    gzFile fp2 = gzopen(r2fq, "r");
    kseq_t *seq1 = kseq_init(fp1);
    kseq_t *seq2 = kseq_init(fp2);
    int l1, l2;
    uint64_t bin;
    int count = 0;
    char pass_fail;
    int blen1_2 = settings->blen / 2;
    int n_len = blen1_2 + settings->homing_sequence_length + settings->offset;
    char *barcode = (char *)malloc((settings->blen + 1) * sizeof(char));
    barcode[settings->blen] = '\0'; // Null-terminate
    l1 = kseq_read(seq1);
    l2 = kseq_read(seq2);
    if(l1 < 0 || l2 < 0) {
            fprintf(stderr, "Could not open fastqs for reading. Abort!\n");
            FREE_MSSI_SETTINGS_PTR(settings);
            FREE_SPLITTER_PTR(splitter);
            exit(EXIT_FAILURE);
    }
    tmp_mseq_t tmp = init_tmp_mseq(seq1->seq.l, settings->blen);
    // Get first barcode.
    set_barcode(seq1, seq2, barcode, settings->offset, blen1_2);
    pass_fail = test_homing_seq(seq1, seq2, settings) ? test_hp_inline(barcode, settings->blen, settings->hp_threshold) : '0';
    mseq_t mvar1 = init_rescale_revcmp_mseq(seq1, barcode, settings->rescaler, &tmp, n_len, 0);
    mseq_t mvar2 = init_rescale_revcmp_mseq(seq2, barcode, settings->rescaler, &tmp, n_len, 1);
    bin = get_binnerul(barcode, settings->n_nucs);
    mseq2fq_inline(splitter->tmp_out_handles_r1[bin], &mvar1, pass_fail);
    mseq2fq_inline(splitter->tmp_out_handles_r2[bin], &mvar2, pass_fail);
    do {
        if(++count % settings->notification_interval == 0) {
            fprintf(stderr, "Number of records processed: %i.\n", count);
        }
        // Iterate through second fastq file.
        set_barcode(seq1, seq2, barcode, settings->offset, blen1_2);
        pass_fail = test_homing_seq(seq1, seq2, settings) ? test_hp_inline(barcode, settings->blen, settings->hp_threshold) : '0';
        //fprintf(stdout, "Randomly testing to see if the reading is working. %s", seq1->seq.s);
        update_mseq(&mvar1, barcode, seq1, settings->rescaler, &tmp, n_len, 0);
        update_mseq(&mvar2, barcode, seq2, settings->rescaler, &tmp, n_len, 1);
        bin = get_binnerul(barcode, settings->n_nucs);
        mseq2fq_inline(splitter->tmp_out_handles_r1[bin], &mvar1, pass_fail);
        mseq2fq_inline(splitter->tmp_out_handles_r2[bin], &mvar2, pass_fail);
    } while (((l1 = kseq_read(seq1)) >= 0) && ((l2 = kseq_read(seq2)) >= 0));
    for(int i = 0; i < splitter->n_handles; ++i) {
        fclose(splitter->tmp_out_handles_r1[i]);
        fclose(splitter->tmp_out_handles_r2[i]);
    }
    tmp_mseq_destroy(tmp);
    mseq_destroy(&mvar1);
    mseq_destroy(&mvar2);
    free(barcode);
    kseq_destroy(seq1);
    kseq_destroy(seq2);
    gzclose(fp1);
    gzclose(fp2);
    fprintf(stderr, "Successfully completed pp_split_inline.\n");
    return splitter;
}

// Rescaling
// TODO:
// 1. Write modified splitmark_core_inline
// 2. Add parser for recalibrated quality scores.


int main(int argc, char *argv[])
{
    // Build settings struct
    int hp_threshold;
    int n_nucs;
    char *output_basename;
    const char *default_basename = "metasyntactic_var";
    int blens[10];
    mssi_settings_t settings = {
        .hp_threshold = 10,
        .n_nucs = 4,
        .output_basename = NULL,
        .input_r1_path = NULL,
        .input_r2_path = NULL,
        .homing_sequence = NULL,
        .n_handles = 0,
        .notification_interval = 100000,
        .blen = 0,
        .homing_sequence_length = 0,
        .offset = 0,
        .rescaler = NULL,
        .rescaler_path = NULL,
        .run_hash_dmp = NULL,
        .ffq_prefix = NULL,
        .threads = 1
    };
    omp_set_dynamic(0); // Tell omp that I want to set my number of threads 4realz
    int c;
    while ((c = getopt(argc, argv, "t:ho:n:s:l:m:r:dp:f:")) > -1) {
        switch(c) {
            case 'n': settings.n_nucs = atoi(optarg); break;
            case 't': settings.hp_threshold = atoi(optarg); break;
            case 'o': settings.output_basename = strdup(optarg); break;
            case 'r': settings.rescaler_path = strdup(optarg); break;
            case 's': settings.homing_sequence = strdup(optarg); settings.homing_sequence_length = strlen(settings.homing_sequence); break;
            case 'l': settings.blen = 2 * atoi(optarg); break;
            case 'm': settings.offset = atoi(optarg); break;
            case 'p': settings.threads = atoi(optarg); omp_set_num_threads(settings.threads); break;
            case 'd': settings.run_hash_dmp = 1; break;
            case 'f': settings.ffq_prefix = strdup(optarg); break;
            case 'h': print_crms_usage(argv); return 0;
            default: print_crms_opt_err(argv, optarg);
        }
    }
    settings.n_handles = ipow(4, settings.n_nucs);
    if(settings.n_handles > get_fileno_limit()) {
        increase_nofile_limit(kroundup32(settings.n_handles));
    }
    fprintf(stderr, "Starting main.\n");
    if(argc < 5) {
        print_crms_usage(argv); exit(1);
    }

    if(settings.ffq_prefix && !settings.run_hash_dmp) {
        fprintf(stderr, "Final fastq prefix option provided but run_hash_dmp not selected."
                "Either eliminate the -f flag or add the -d flag.\n");
        exit(EXIT_FAILURE);
    }

    if(!settings.homing_sequence) {
        fprintf(stderr, "Homing sequence not provided. Required.\n");
        exit(EXIT_FAILURE);
    }
    if(!settings.blen) {
        fprintf(stderr, "Barcode length not provided. Required. Abort!\n");
        exit(EXIT_FAILURE);
    }

    if(settings.offset) {
        settings.blen -= 2 * settings.offset;
    }
    if(settings.rescaler_path) {
        settings.rescaler = parse_rescaler(settings.rescaler_path);
    }
    fprintf(stderr, "About to get the read paths.\n");
    char r1fq[100];
    char r2fq[100];
    if(argc - 1 != optind + 1) {
        fprintf(stderr, "Both read 1 and read 2 fastqs are required. See usage.\n", argc, optind);
        print_crms_usage(argv);
        return 1;
    }
    strcpy(r1fq, argv[optind]);
    strcpy(r2fq, argv[optind + 1]);

    if(!settings.output_basename) {
        settings.output_basename = make_crms_outfname(r1fq);
        fprintf(stderr, "Output basename not provided. Defaulting to variation on input: %s.\n", settings.output_basename);
    }
    mark_splitter_t *splitter = pp_split_inline1(r1fq, r2fq, &settings);
    if(settings.rescaler) {
        cfree_rescaler(settings);
    }
    if(settings.run_hash_dmp) {
        fprintf(stderr, "Now executing hash dmp.\n");
        if(!settings.ffq_prefix) {
            settings.ffq_prefix = make_default_outfname(r1fq, ".dmp.final");
        }
        // Whatever I end up putting into here.
        splitterhash_params_t *params = init_splitterhash(&settings, splitter);
        #pragma omp parallel
        {
            #pragma omp for
            for(int i = 0; i < settings.n_handles; ++i) {
                fprintf(stderr, "Now running omgz core on input filename %s and output filename %s.\n",
                        params->infnames_r1[i], params->outfnames_r1[i]);
                omgz_core(params->infnames_r1[i], params->outfnames_r1[i]);
                omgz_core(params->infnames_r2[i], params->outfnames_r2[i]);
            }
        }
        // Remove temporary split files
        fprintf(stderr, "Now removing temporary files.\n");
        char del_buf[250];
        for(int i = 0; i < splitter->n_handles; ++i) {
            fprintf(stderr, "Now removing temporary files %s and %s.\n",
                    splitter->fnames_r1[i], splitter->fnames_r2[i]);
            sprintf(del_buf, "rm %s %s", splitter->fnames_r1[i], splitter->fnames_r2[i]);
            system(del_buf);
        }
        fprintf(stderr, "Now building cat string.\n");
        char cat_buff1[CAT_BUFFER_SIZE] = "/bin/cat ";
        char cat_buff2[CAT_BUFFER_SIZE] = "/bin/cat ";
        fprintf(stderr, "Current cat buff1: %s.\n", cat_buff1);
        for(int i = 0; i < settings.n_handles; ++i) {
            strcat(cat_buff1, params->outfnames_r1[i]);
            strcat(cat_buff1, " ");
            strcat(cat_buff2, params->outfnames_r2[i]);
            strcat(cat_buff2, " ");
        }
        strcat(cat_buff1, " > ");
        strcat(cat_buff1, settings.ffq_prefix);
        strcat(cat_buff1, ".R1.fq");
        strcat(cat_buff2, " > ");
        strcat(cat_buff2, settings.ffq_prefix);
        strcat(cat_buff2, ".R2.fq");
        fprintf(stderr, "Now calling cat string '%s'.\n", cat_buff1);
        system(cat_buff1);
        fprintf(stderr, "Now calling cat string '%s'.\n", cat_buff2);
        system(cat_buff2);
        for(int i = 0; i < splitter->n_handles; ++i) {
            fprintf(stderr, "Now calling 'rm %s %s'\n", params->outfnames_r1[i], params->outfnames_r2[i]);
            sprintf(del_buf, "rm %s %s", params->outfnames_r1[i], params->outfnames_r2[i]);
            system(del_buf);
        }
        splitterhash_destroy(params);
        free(settings.ffq_prefix);
    }
    free_mssi_settings(settings);
    FREE_SPLITTER_PTR(splitter);
    return 0;
}
