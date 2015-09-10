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
#include "fqmarksplit.h"

// Inline function declarations
//int lh3_sort_call(char *fname, char *outfname);
void FREE_SPLITTER(mark_splitter_t var);

//void apply_lh3_sorts(sort_overlord_t *dispatcher, mss_settings_t *settings);
int ipow(int base, int exp);
mark_splitter_t init_splitter(mss_settings_t *settings_ptr);
int get_binner(char binner[], int length);
char test_hp_inline(char *barcode, int length, int threshold);
char test_hp(kseq_t *seq, int threshold);

// Macros



void print_usage(char *argv[])
{
        fprintf(stderr, "Usage: %s <options> <Fq.R1.seq> <Fq.R2.seq>"
                        "\nFlags:\n"
                        "-h: Print usage.\n"
                        "-t: Homopolymer failure threshold. A molecular barcode with"
                        " a homopolymer of length >= this limit is flagged as QC fail."
                        "Default: 10.\n"
                        "-o: Output basename. Currently required, as string "
                        "manipulation in C is a bit of work and I'd rather spend my "
                        "time building code than messing around with string "
                        "manipulation that doesn't add to the code base.\n"
                        "-n: Number of nucleotides at the beginning of the barcode to use to split the output.\n"
                        "-m: Mask first n nucleotides in read for barcode. Default: 0. Recommended: 1.\n"
                        "-l: Number of nucleotides at the beginning of each read to "
                        "use for barcode. Final barcode length is twice this.\n"
                        "-s: homing sequence. If not provided, %s will not look for it.\n", argv[0], argv[0]);
}

void print_opt_err(char *argv[], char *optarg)
{
    print_usage(argv);
    fprintf(stderr, "Unrecognized option %s. Abort!\n", optarg);
    exit(1);
}


int test_homing_seq(kseq_t *seq1, kseq_t *seq2, mssi_settings_t *settings_ptr) {
    if(!settings_ptr->homing_sequence) {
        return 1;
    }
    else {
        return (memcmp(seq1 + (settings_ptr->blen / 2),
                       settings_ptr->homing_sequence,
                       settings_ptr->homing_sequence_length) == 0 &&
                memcmp(seq2 + (settings_ptr->blen / 2),
                       settings_ptr->homing_sequence,
                       settings_ptr->homing_sequence_length) == 0);
    }
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

void splitmark_core_inline(kseq_t *seq1, kseq_t *seq2,
                           mssi_settings_t settings, mark_splitter_t splitter)
{
    int l1, l2, bin;
    int count = 0;
    char pass_fail;
    int readlen = 0;
    int n_len = settings.blen + settings.homing_sequence_length;
    char * barcode;
    char * tmp_n_str;
    int blen1_2 = settings.blen / 2;
    barcode = (char *)malloc((settings.blen + 1) * sizeof(char));
    barcode[settings.blen] = '\0'; // Null-terminate
    l1 = kseq_read(seq1);
    if(l1 < 0) {
            fprintf(stderr, "Could not open fastq for reading. Abort!\n");
            FREE_MSSI_SETTINGS(settings);
            FREE_SPLITTER(splitter);
            exit(EXIT_FAILURE);
    }
    readlen = strlen(seq1->seq.s);
    tmp_n_str = (char *)malloc((readlen + 1) * sizeof(char));
    tmp_n_str[readlen] = '\0';
    do {
        count += 1;
        if(!(count % settings.notification_interval)) {
            fprintf(stderr, "Number of records processed: %i.\n", count);
        }
        // Iterate through second fastq file.
        l2 = kseq_read(seq2);
        if (l2 < 0) {
            fprintf(stderr, "Read 2 return value for kseq_read less than "
                            "0. Are the fastqs different sizes? Abort!\n");
            FREE_MSSI_SETTINGS(settings);
            FREE_SPLITTER(splitter);
            exit(EXIT_FAILURE);
        }
        memcpy(barcode, seq1->seq.s + settings.offset, blen1_2 * sizeof(char)); // Copying the fist half of the barcode
        memcpy(barcode + blen1_2, seq2->seq.s + settings.offset,
               blen1_2 * sizeof(char));
        pass_fail = test_homing_seq(seq1, seq2, &settings) ? test_hp_inline(barcode, settings.blen, settings.hp_threshold) : '0';
        //fprintf(stdout, "Randomly testing to see if the reading is working. %s", seq1->seq.s);
        bin = get_binner(barcode, settings.n_nucs);
        KSEQ_2_FQ_INLINE(splitter.tmp_out_handles_r1[bin], seq1, barcode, pass_fail, tmp_n_str, readlen, n_len);
        KSEQ_2_FQ_INLINE(splitter.tmp_out_handles_r2[bin], seq2, barcode, pass_fail, tmp_n_str, readlen, n_len);
    } while ((l1 = kseq_read(seq1)) >= 0);
    free(barcode);
    free(tmp_n_str);
}

// Rescaling
// TODO:
// 1. Write modified splitmark_core_inline
// 2. Add parser for recalibrated quality scores.


void splitmarkrecalibratecrc_inline(kseq_t *seq1, kseq_t *seq2,
                                    mssi_settings_t settings, mark_splitter_t splitter)
{
    int l1, l2, bin;
    int count = 0;
    char pass_fail;
    int readlen = 0;
    int n_len = settings.blen + settings.homing_sequence_length;
    char * barcode;
    char * tmp_n_str;
    int blen1_2 = settings.blen / 2;
    barcode = (char *)malloc((settings.blen + 1) * sizeof(char));
    barcode[settings.blen] = '\0'; // Null-terminate
    l1 = kseq_read(seq1);
    if(l1 < 0) {
            fprintf(stderr, "Could not open fastq for reading. Abort!\n");
            FREE_MSSI_SETTINGS(settings);
            FREE_SPLITTER(splitter);
            exit(EXIT_FAILURE);
    }
    readlen = strlen(seq1->seq.s);
    tmp_n_str = (char *)malloc((readlen + 1) * sizeof(char));
    tmp_n_str[readlen] = '\0';
    do {
        count += 1;
        if(!(count % settings.notification_interval)) {
            fprintf(stderr, "Number of records processed: %i.\n", count);
        }
        // Iterate through second fastq file.
        l2 = kseq_read(seq2);
        if (l2 < 0) {
            fprintf(stderr, "Read 2 return value for kseq_read less than "
                            "0. Are the fastqs different sizes? Abort!\n");
            FREE_MSSI_SETTINGS(settings);
            FREE_SPLITTER(splitter);
            exit(EXIT_FAILURE);
        }
        memcpy(barcode, seq1->seq.s + settings.offset, blen1_2 * sizeof(char)); // Copying the fist half of the barcode
        memcpy(barcode + blen1_2, seq2->seq.s + settings.offset,
               blen1_2 * sizeof(char));
        pass_fail = test_homing_seq(seq1, seq2, &settings) ? test_hp_inline(barcode, settings.blen, settings.hp_threshold) : '0';
        //fprintf(stdout, "Randomly testing to see if the reading is working. %s", seq1->seq.s);
        bin = get_binner(barcode, settings.n_nucs);
        KSEQ_2_FQ_INLINE(splitter.tmp_out_handles_r1[bin], seq1, barcode, pass_fail, tmp_n_str, readlen, n_len);
        KSEQ_2_FQ_INLINE(splitter.tmp_out_handles_r2[bin], seq2, barcode, pass_fail, tmp_n_str, readlen, n_len);
    } while ((l1 = kseq_read(seq1)) >= 0);
    free(barcode);
    free(tmp_n_str);
}


int main(int argc, char *argv[])
{
    // Build settings struct
    int hp_threshold;
    int n_nucs;
    char *output_basename;
    int threads;
    const char *default_basename = "metasyntactic_var";
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
        .rescaler_path = NULL
    };
    int c;
    while ((c = getopt(argc, argv, "t:ho:n:s:l:m:r:")) > -1) {
        switch(c) {
            case 'n': settings.n_nucs = atoi(optarg); break;
            case 't': settings.hp_threshold = atoi(optarg); break;
            case 'o': settings.output_basename = strdup(optarg); break;
            case 'r': settings.rescaler_path = strdup(optarg); break;
            case 's': settings.homing_sequence = strdup(optarg); settings.homing_sequence_length = strlen(settings.homing_sequence); break;
            case 'l': settings.blen = 2 * atoi(optarg); break;
            case 'm': settings.offset = atoi(optarg); break;
            case 'h': print_usage(argv); return 0;
            default: print_opt_err(argv, optarg);
        }
    }
    settings.n_handles = ipow(4, settings.n_nucs);

    if(!settings.homing_sequence) {
        fprintf(stderr, "Homing sequence not provided. Will not check for it. "
                        "Reads will not be QC failed for its absence.\n");
    }
    if(!settings.blen) {
        fprintf(stderr, "Barcode length not provided. Required. Abort!\n.");
        exit(EXIT_FAILURE);
    }

    if(settings.offset) {
        settings.blen -= 2 * settings.offset;
    }
    if(settings.rescaler_path) {
        settings.rescaler = parse_rescaler(settings.rescaler_path);
    }

    if(!settings.output_basename) {
        fprintf(stderr, "Output basename required. See usage.\n");
        print_usage(argv);
        return 1;
    }
    fprintf(stderr, "About to get the read paths.\n");
    char r1_fq_buf[100];
    char r2_fq_buf[100];
    if(argc - 1 != optind + 1) {
        fprintf(stderr, "Both read 1 and read 2 fastqs are required. See usage.\n", argc, optind);
        print_usage(argv);
        return 1;
    }
    strcpy(r1_fq_buf, argv[optind]);
    strcpy(r2_fq_buf, argv[optind + 1]);
    mssi_settings_t *settings_ptr = &settings;
    gzFile fp_read1, fp_read2;
    fp_read1 = gzopen(r1_fq_buf, "r");
    fp_read2 = gzopen(r2_fq_buf, "r");
    kseq_t *seq1;
    kseq_t *seq2;
    int l1, l2, l_index;
    seq1 = kseq_init(fp_read1);
    seq2 = kseq_init(fp_read2);
    mark_splitter_t splitter = init_splitter_inline(settings_ptr);
    mark_splitter_t *splitter_ptr = &splitter;
    splitmark_core_inline(seq1, seq2,
                          settings, splitter);
    free_mssi_settings(settings);
    FREE_SPLITTER(splitter);
    kseq_destroy(seq1);
    kseq_destroy(seq2);
    return 0;
}
