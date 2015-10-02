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



inline blens_t *get_blens(char *str2parse)
{
    blens_t *ret = (blens_t *)malloc(sizeof(blens_t));
    ret->min_blen = 200;
    ret->max_blen = 0;
    ret->n = 0;
    ret->current_blen = -1;
    char *token;
    while((token = strtok(str2parse, ",")) != NULL){
       ret->blens[ret->n] = atoi(token);
       if(ret->max_blen < ret->blens[ret->n]) {
           ret->max_blen = ret->blens[ret->n];
       }
       else if(ret->min_blen > ret->blens[ret->n]) {
           ret->min_blen = ret->blens[ret->n];
       }
       ++ret->n;
    }
    return ret;
}


inline void free_crms_settings(crms_settings_t settings)
{
    free(settings.output_basename);
    free(settings.input_r1_path);
    free(settings.input_r2_path);
    cond_free(settings.rescaler_path);
    // Note: rescaler array should be handled elsewhere!
    return;
}

void print_usage(char *argv[])
{
        fprintf(stderr, "Usage: %s <options> <Fq.R1.seq> <Fq.R2.seq>"
                        "\nFlags:\n"
                        "-s: homing sequence. REQUIRED.\n"
                        "-l: Number of nucleotides at the beginning of each read to "
                        "use for barcode. Final barcode length is twice this. REQUIRED.\n"
                        "-o: Output basename. Defaults to a variation on input filename.\n"
                        "-t: Homopolymer failure threshold. A molecular barcode with"
                        " a homopolymer of length >= this limit is flagged as QC fail."
                        "Default: 10.\n"
                        "-n: Number of nucleotides at the beginning of the barcode to use to split the output. Default: 4.\n"
                        "-m: Mask first n nucleotides in read for barcode. Default: 0. Recommended: 1.\n"
                        "-h: Print usage.\n", argv[0], argv[0]);
}

void print_opt_err(char *argv[], char *optarg)
{
    print_usage(argv);
    fprintf(stderr, "Unrecognized option %s. Abort!\n", optarg);
    exit(1);
}


mark_splitter_t init_splitter_crms(crms_settings_t* settings_ptr)
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

inline int vl_homing_loc(kseq_t *seq1, kseq_t *seq2, crms_settings_t *settings_ptr)
{
#if !NDEBUG
    if(!settings_ptr->blen_data->homing_sequence) {
        fprintf(stderr, "Okay, if the homing sequence isn't set, just give up.\n");
        exit(EXIT_FAILURE);
    }
#endif
    for(int i = 0; i < settings_ptr->blen_data->n; ++i) {
        if(memcmp(seq1->seq.s + settings_ptr->blen_data->blens[i] + settings_ptr->offset,
                  settings_ptr->blen_data->homing_sequence,
                  settings_ptr->blen_data->homing_sequence_length) == 0) {
            return settings_ptr->blen_data->blens[i];
        }
    }
    return 0;
}


/*
 * Pre-processes and splits fastqs with inline barcodes.
 */
void vl_split_inline(kseq_t *seq1, kseq_t *seq2,
                     crms_settings_t *settings, mark_splitter_t *splitter)
{
    int l1, l2;
    uint64_t bin;
#if LOGGING
    int count = 0;
#endif
    char pass_fail;
    int readlen = 0;
    char *barcode;
    barcode = (char *)malloc((settings->blen_data->max_blen + 1) * sizeof(char));
    barcode[settings->blen_data->max_blen] = '\0'; // Null-terminate
    l1 = kseq_read(seq1);
    l2 = kseq_read(seq2);
    if(l1 < 0 || l2 < 0) {
            fprintf(stderr, "Could not open fastqs for reading. Abort!\n");
            FREE_SPLITTER(*splitter);
            exit(EXIT_FAILURE);
    }
    readlen = seq1->seq.l;
    settings->blen_data->current_blen =  vl_homing_loc(seq1, seq2, settings);
    int n_len = settings->blen_data->current_blen + settings->blen_data->homing_sequence_length + settings->offset;
    tmp_mseq_t tmp = init_tmp_mseq(readlen, settings->blen_data->max_blen);
    // Get first barcode.
    set_barcode(seq1, seq2, barcode, settings->offset, settings->blen_data->min_blen);
    pass_fail = settings->blen_data->current_blen ? test_hp_inline(barcode, settings->blen_data->current_blen, settings->hp_threshold) : '0';
    mseq_t mvar1 = init_rescale_revcmp_mseq(seq1, barcode, settings->rescaler, &tmp, n_len, 0);
    mseq_t mvar2 = init_rescale_revcmp_mseq(seq2, barcode, settings->rescaler, &tmp, n_len, 1);
    bin = get_binnerul(barcode, settings->n_nucs);
    mseq2fq_inline(splitter->tmp_out_handles_r1[bin], &mvar1, pass_fail);
    mseq2fq_inline(splitter->tmp_out_handles_r2[bin], &mvar2, pass_fail);
    do {
        if(!(count++ % settings->notification_interval)) fprintf(stderr, "Number of records processed: %i.\n", count);
        // Iterate through second fastq file.
        set_barcode(seq1, seq2, barcode, settings->offset, settings->blen_data->min_blen);
        settings->blen_data->current_blen = vl_homing_loc(seq1, seq2, settings);
        pass_fail = settings->blen_data->current_blen ? test_hp_inline(barcode, settings->blen_data->current_blen, settings->hp_threshold) : '0';
        n_len = settings->blen_data->current_blen + settings->blen_data->homing_sequence_length + settings->offset;
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
    return;
}


int main(int argc, char *argv[])
{
    if(argc < 5) {
        print_usage(argv); exit(1);
    }
    // Build settings struct
    int hp_threshold;
    int n_nucs;
    char *output_basename;
    int threads;
    const char *default_basename = "metasyntactic_var";
    crms_settings_t settings = {
        .hp_threshold = 10,
        .n_nucs = 4,
        .output_basename = NULL,
        .input_r1_path = NULL,
        .input_r2_path = NULL,
        .n_handles = 0,
        .notification_interval = 1000000,
        .blen_data = (blens_t *)calloc(1, sizeof(blens_t)),
        .offset = 0,
        .rescaler = NULL,
        .rescaler_path = NULL
    };
    settings.blen_data->max_blen = -1;
    settings.blen_data->homing_sequence_length = 0;
    int c;
    while ((c = getopt(argc, argv, "t:ho:n:s:l:m:r:")) > -1) {
        switch(c) {
            case 'n': settings.n_nucs = atoi(optarg); break;
            case 't': settings.hp_threshold = atoi(optarg); break;
            case 'o': settings.output_basename = strdup(optarg); break;
            case 'r': settings.rescaler_path = strdup(optarg); break;
            case 's': strcpy(settings.blen_data->homing_sequence, optarg); settings.blen_data->homing_sequence_length = strlen(settings.blen_data->homing_sequence); break;
            case 'l': settings.blen_data = get_blens(optarg); break;
            case 'm': settings.offset = atoi(optarg); break;
            case 'h': print_usage(argv); return 0;
            default: print_opt_err(argv, optarg);
        }
    }
    if(!settings.blen_data->homing_sequence_length) {
        fprintf(stderr, "homing sequence not provided. e.g., %s -s <homing_sequence> <other_args>.\n", argv[0]);
        exit(EXIT_FAILURE);

    }
    if(settings.blen_data->max_blen < 0) {
        fprintf(stderr, "blen data not provided. This should be a comma-separated list "
                         "of positive integers for possible barcode lengths.\n");
        fprintf(stderr, "e.g., %s -l 8,9,10 <other_args>.\n", argv[0]);
        exit(EXIT_FAILURE);
    }
    settings.n_handles = ipow(4, settings.n_nucs);
    if(settings.n_handles > get_fileno_limit()) {
        increase_nofile_limit(kroundup32(settings.n_handles));
    }
    fprintf(stderr, "Starting main.\n");

    if(settings.offset) {
        for(int i = 0; i < settings.blen_data->n; ++i) {
            settings.blen_data->blens[i] -= 2 * settings.offset;
        }
        settings.blen_data->max_blen -= 2 * settings.offset;
        settings.blen_data->min_blen -= 2 * settings.offset;
    }
    if(settings.rescaler_path) {
        settings.rescaler = parse_rescaler(settings.rescaler_path);
    }
    fprintf(stderr, "About to get the read paths.\n");
    char r1fq[100];
    char r2fq[100];
    if(argc - 1 != optind + 1) {
        fprintf(stderr, "Both read 1 and read 2 fastqs are required. See usage.\n", argc, optind);
        print_usage(argv);
        return 1;
    }
    strcpy(r1fq, argv[optind]);
    strcpy(r2fq, argv[optind + 1]);

    if(!settings.output_basename) {
        fprintf(stderr, "Output basename not provided. Defaulting to variation on input: %s.\n", settings.output_basename);
        settings.output_basename = make_crms_outfname(r1fq);
    }
    gzFile fp_read1, fp_read2;
    fp_read1 = gzopen(r1fq, "r");
    fp_read2 = gzopen(r2fq, "r");
    int l1, l2;
    kseq_t *seq1 = kseq_init(fp_read1);
    kseq_t *seq2 = kseq_init(fp_read2);
    mark_splitter_t splitter = init_splitter_crms(&settings);
    vl_split_inline(seq1, seq2, &settings, &splitter);
    cfree_rescaler(settings);
    free(settings.blen_data);
    free_crms_settings(settings);
    FREE_SPLITTER(splitter);
    kseq_destroy(seq1);
    kseq_destroy(seq2);
    return 0;
}
