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
#include "vlcrms.h"

#define MAX_N_BLENS 5
#define MAX_HOMING_SEQUENCE 8

typedef struct blens {
	int max_blen; // Last value in blens
	int min_blen; // Lowest value in blens
	int blens[MAX_N_BLENS]; // Array holding blens
	int n; // Number of blens to look for
	int current_blen;
	int homing_sequence_length;
	char homing_sequence[MAX_HOMING_SEQUENCE + 1];
} blens_t;

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


inline int get_nlen()

inline void free_crms_settings(crms_settings_t settings)
{
    free(settings.output_basename);
    free(settings.input_r1_path);
    free(settings.input_r2_path);
    if(settings.rescaler_path) free(settings.rescaler_path);
    // Note: rescaler array should be handled elsewhere!
    return;
}

typedef struct crms_settings {
    int hp_threshold; // The minimum length of a homopolymer run to fail a barcode.
    int n_nucs; // Number of nucleotides to split by.
    char *output_basename;
    char *input_r1_path;
    char *input_r2_path;
    char *homing_sequence; // Homing sequence...
    int homing_sequence_length; // Length of homing sequence, should it be used.
    int n_handles; // Number of handles
    int notification_interval; // How many sets of records do you want to process between progress reports?
    blens_t *blen_data;
    int offset; // Number of bases at the start of the inline barcodes to skip for low quality.
    char ****rescaler; // Three-dimensional rescaler array. Size: [readlen, 39, 4] (length of reads, number of original quality scores, number of bases)
    char *rescaler_path; // Path to rescaler for
} crms_settings_t;

// Inline function declarations
int crc_flip(mseq_t *mvar, char *barcode, int blen, int readlen);
void crc_mseq(mseq_t *mvar, tmp_mseq_t *tmp);
void FREE_SPLITTER(mark_splitter_t var);
mseq_t init_rescale_revcmp_mseq(kseq_t *seq, char *barcode, char ****rescaler, tmp_mseq_t *tmp, int n_len, int is_read2);
mark_splitter_t init_splitter_inline(mssi_settings_t* settings_ptr);
mark_splitter_t init_splitter(mss_settings_t *settings_ptr);
int ipow(int base, int exp);
void mseq_destroy(mseq_t *mvar);
void mseq_rescale_init(kseq_t *seq, mseq_t *ret, char ****rescaler, tmp_mseq_t *tmp, int n_len, int is_read2);
int nuc2num(char character);
int nuc_cmp(char forward, char reverse);
char ****parse_rescaler(char *qual_rescale_fname);
int rescale_qscore(int readnum, int qscore, int cycle, char base, char ****rescaler);
void set_barcode(kseq_t *seq1, kseq_t *seq2, char *barcode, int offset, int blen1_2);
int test_homing_seq(kseq_t *seq1, kseq_t *seq2, mssi_settings_t *settings_ptr);
char test_hp_inline(char *barcode, int length, int threshold);
char test_hp(kseq_t *seq, int threshold);
void tmp_mseq_destroy(tmp_mseq_t mvar);
void update_mseq(mseq_t *mvar, char *barcode, kseq_t *seq, char ****rescaler, tmp_mseq_t *tmp, int n_len, int is_read2);
char nuc_cmpl(char character);
void mseq2fq_inline(FILE *handle, mseq_t *mvar, char pass_fail);
int count_lines(char *fname);
void FREE_SPLITTER(mark_splitter_t var);
char *trim_ext(char *fname);
char *make_default_outfname(char *fname, const char *suffix);
char *mark_crms_outfname(char *fname);
void u32toa_branchlut(uint32_t value, char* buffer);
void i32toa_branchlut(int32_t value, char* buffer);
int get_fileno_limit();
void increase_nofile_limit(int new_limit);
uint64_t get_binnerul(char *barcode, int length);
int get_binner(char *barcode, int length);
uint64_t ulpow(uint64_t base, uint64_t exp);


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

/*
 * Pre-processes and splits fastqs with inline barcodes.
 */
void vl_split_inline(kseq_t *seq1, kseq_t *seq2,
                     crms_settings_t settings, mark_splitter_t splitter)
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
        if(!(count++ % settings.notification_interval)) fprintf(stderr, "Number of records processed: %i.\n", count);
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
    tmp_mseq_destroy(tmp);
    mseq_destroy(&mvar1);
    mseq_destroy(&mvar2);
    free(barcode);
    return;
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
    int threads;
    const char *default_basename = "metasyntactic_var";
    crms_settings_t settings = {
        .hp_threshold = 10,
        .n_nucs = 4,
        .output_basename = NULL,
        .input_r1_path = NULL,
        .input_r2_path = NULL,
        .homing_sequence = NULL,
        .n_handles = 0,
        .notification_interval = 100000,
        .blen_data = (blens_t *)malloc(sizeof(blens_t)),
        .homing_sequence_length = 0,
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
            case 's': settings.blen_data->homing_sequence = strdup(optarg); settings.blen_data->homing_sequence_length = strlen(settings.homing_sequence); break;
            case 'l': settings.blen_data = get_blens(optarg); break;
            case 'm': settings.offset = atoi(optarg); break;
            case 'h': print_usage(argv); return 0;
            default: print_opt_err(argv, optarg);
        }
    }
    if(!settings.blen_data) {
    	settings.blen_data->homing_sequence_length
    }
    if(!settings.blen_data->homing_sequence_length) {
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
    if(argc < 5) {
        print_usage(argv); exit(1);
    }

    if(!settings.homing_sequence) {
        fprintf(stderr, "Homing sequence not provided. Will not check for it. "
                        "Reads will not be QC failed for its absence.\n");
        fprintf(stderr, "Actually, that's stupid. I'm going to exit. Abort!\n");
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
    pp_split_inline(seq1, seq2, settings, splitter);
    if(settings.rescaler) {
        int readlen = count_lines(settings.rescaler_path);
        for(int i = 0; i < 2; ++i) {
            for(int j = 0; j < readlen; ++j) {
                for(int k = 0; k < 39; ++k) {
                	if(settings.rescaler[i][j][k]) {
                       free(settings.rescaler[i][j][k]);
                       settings.rescaler[i][j][k] = NULL;
                	}
                }
                if(settings.rescaler[i][j]) {
                    free(settings.rescaler[i][j]);
                    settings.rescaler[i][j] = NULL;
                }
            }
            if(settings.rescaler[i]) {
                free(settings.rescaler[i]);
                settings.rescaler[i] = NULL;
            }
        }
        free(settings.rescaler);
        settings.rescaler = NULL;
    }
    free_crms_settings(settings);
    FREE_SPLITTER(splitter);
    kseq_destroy(seq1);
    kseq_destroy(seq2);
    return 0;
}
.blen_data->
