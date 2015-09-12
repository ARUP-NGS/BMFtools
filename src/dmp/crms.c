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

// Inline function declarations
int crc_flip(mseq_t *mvar, char *barcode, int blen, int readlen);
void crc_mseq(mseq_t *mvar, tmp_mseq_t *tmp);
void FREE_SPLITTER(mark_splitter_t var);
int get_binner(char binner[], int length);
mseq_t init_rescale_revcmp_mseq(kseq_t *seq, char *barcode, char ***rescaler, tmp_mseq_t *tmp, int n_len);
mark_splitter_t init_splitter_inline(mssi_settings_t* settings_ptr);
mark_splitter_t init_splitter(mss_settings_t *settings_ptr);
int ipow(int base, int exp);
void mseq_destroy(mseq_t *mvar);
void mseq_rescale_init(kseq_t *seq, mseq_t *ret, char ***rescaler, tmp_mseq_t *tmp, int n_len);
int nuc2num(char character);
int nuc_cmp(char forward, char reverse);
char ***parse_rescaler(char *qual_rescale_fname);
int rescale_qscore(int qscore, int cycle, char base, char ***rescaler);
void set_barcode(kseq_t *seq1, kseq_t *seq2, char *barcode, int offset, int blen1_2);
int test_homing_seq(kseq_t *seq1, kseq_t *seq2, mssi_settings_t *settings_ptr);
char test_hp_inline(char *barcode, int length, int threshold);
char test_hp(kseq_t *seq, int threshold);
void tmp_mseq_destroy(tmp_mseq_t mvar);
void update_mseq(mseq_t *mvar, char *barcode, kseq_t *seq, char ***rescaler, tmp_mseq_t *tmp, int n_len);
char nuc_cmpl(char character);
void mseq2fq_inline(FILE *handle, mseq_t *mvar, char pass_fail);
int count_lines(char *fname);
void FREE_SPLITTER(mark_splitter_t var);
char *** parse_rescaler(char *qual_rescale_fname);
char *trim_ext(char *fname);
char *make_default_outfname(char *fname, const char *suffix);
char *mark_crms_outfname(char *fname);


void print_usage(char *argv[])
{
        fprintf(stderr, "Usage: %s <options> <Fq.R1.seq> <Fq.R2.seq>"
                        "\nFlags:\n"
                        "-l: Number of nucleotides at the beginning of each read to "
                        "use for barcode. Final barcode length is twice this. REQUIRED.\n"
                        "-o: Output basename. Defaults to a variation on input filename.\n"
                        "-t: Homopolymer failure threshold. A molecular barcode with"
                        " a homopolymer of length >= this limit is flagged as QC fail."
                        "Default: 10.\n"
                        "-n: Number of nucleotides at the beginning of the barcode to use to split the output. Default: 4.\n"
                        "-m: Mask first n nucleotides in read for barcode. Default: 0. Recommended: 1.\n"
                        "-s: homing sequence. If not provided, %s will not look for it.\n"
                        "-h: Print usage.\n", argv[0], argv[0]);
}

void print_opt_err(char *argv[], char *optarg)
{
    print_usage(argv);
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
 * Pre-processes and splits fastqs with inline barcodes.
 */
void pp_split_inline(kseq_t *seq1, kseq_t *seq2,
                     mssi_settings_t settings, mark_splitter_t splitter)
{
    int l1, l2, bin;
    int count = 0;
    char pass_fail;
    int readlen = 0;
    char * barcode;
    char * tmp_n_str;
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
    readlen = strlen(seq1->seq.s);
    tmp_n_str = (char *)malloc((readlen + 1) * sizeof(char));
    tmp_n_str[readlen] = '\0';
    tmp_mseq_t tmp = init_tmp_mseq(readlen, settings.blen);
    // Get first barcode.
    set_barcode(seq1, seq2, barcode, settings.offset, blen1_2);
    pass_fail = test_homing_seq(seq1, seq2, &settings) ? test_hp_inline(barcode, settings.blen, settings.hp_threshold) : '0';
    mseq_t mvar1 = init_rescale_revcmp_mseq(seq1, barcode, settings.rescaler, &tmp, n_len);
    mseq_t mvar2 = init_rescale_revcmp_mseq(seq2, barcode, settings.rescaler, &tmp, n_len);
    bin = get_binner(barcode, settings.n_nucs);
    mseq2fq_inline(splitter.tmp_out_handles_r1[bin], &mvar1, pass_fail);
    mseq2fq_inline(splitter.tmp_out_handles_r2[bin], &mvar2, pass_fail);
    do {
        count += 1;
        if(!(count % settings.notification_interval)) fprintf(stderr, "Number of records processed: %i.\n", count);
        // Iterate through second fastq file.
        set_barcode(seq1, seq2, barcode, settings.offset, blen1_2);
        pass_fail = test_homing_seq(seq1, seq2, &settings) ? test_hp_inline(barcode, settings.blen, settings.hp_threshold) : '0';
        //fprintf(stdout, "Randomly testing to see if the reading is working. %s", seq1->seq.s);
        update_mseq(&mvar1, barcode, seq1, settings.rescaler, &tmp, n_len);
        update_mseq(&mvar2, barcode, seq2, settings.rescaler, &tmp, n_len);
        bin = get_binner(barcode, settings.n_nucs);
        mseq2fq_inline(splitter.tmp_out_handles_r1[bin], &mvar1, pass_fail);
        mseq2fq_inline(splitter.tmp_out_handles_r2[bin], &mvar2, pass_fail);
    } while (((l1 = kseq_read(seq1)) >= 0) && ((l2 = kseq_read(seq2)) >= 0));
    tmp_mseq_destroy(tmp);
    mseq_destroy(&mvar1);
    mseq_destroy(&mvar2);
    free(barcode);
    free(tmp_n_str);

}

// Rescaling
// TODO:
// 1. Write modified splitmark_core_inline
// 2. Add parser for recalibrated quality scores.


int main(int argc, char *argv[])
{
    fprintf(stderr, "Starting main.\n");
    if(argc < 5) {
        print_usage(argv); exit(1);
    }
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
    fprintf(stderr, "Parsing args.\n");
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
    char r1_fq_buf[100];
    char r2_fq_buf[100];
    if(argc - 1 != optind + 1) {
        fprintf(stderr, "Both read 1 and read 2 fastqs are required. See usage.\n", argc, optind);
        print_usage(argv);
        return 1;
    }
    strcpy(r1_fq_buf, argv[optind]);
    strcpy(r2_fq_buf, argv[optind + 1]);

    if(!settings.output_basename) {
        settings.output_basename = mark_crms_outfname(r1_fq_buf);
        fprintf(stderr, "Output basename not provided. Defaulting to variation on input: %s.\n", settings.output_basename);
    }
    gzFile fp_read1, fp_read2;
    fp_read1 = gzopen(r1_fq_buf, "r");
    fp_read2 = gzopen(r2_fq_buf, "r");
    int l1, l2;
    kseq_t *seq1 = kseq_init(fp_read1);
    kseq_t *seq2 = kseq_init(fp_read2);
    mark_splitter_t splitter = init_splitter_inline(&settings);
    pp_split_inline(seq1,seq2, settings, splitter);
    if(settings.rescaler) {
        int readlen = count_lines(settings.rescaler_path);
        for(int i = 0; i < readlen; i++) {
            for(int j = 0; j < 39; j++) {
                free(settings.rescaler[i][j]);
            }
            free(settings.rescaler[i]);
        }
    }
    free_mssi_settings(settings);
    FREE_SPLITTER(splitter);
    kseq_destroy(seq1);
    kseq_destroy(seq2);
    return 0;
}
