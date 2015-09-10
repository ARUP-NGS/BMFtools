#include <getopt.h>
#include <math.h>
#include <omp.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include "include/kseq.h"

KSEQ_INIT(gzFile, gzread)


typedef struct mss_settings {
    int hp_threshold;
    int n_nucs;
    int n_handles;
    char *output_basename;
    int threads;
    int notification_interval;
    char *index_fq_path;
    char *input_r1_path;
    char *input_r2_path;
} mss_settings_t;


typedef struct mssi_settings {
    int hp_threshold; // The minimum length of a homopolymer run to fail a barcode.
    int n_nucs; // Number of nucleotides to split by.
    char *output_basename;
    char *input_r1_path;
    char *input_r2_path;
    char *homing_sequence; // Homing sequence...
    int homing_sequence_length; // Length of homing sequence, should it be used.
    int n_handles; // Number of handles
    int notification_interval; // How many sets of records do you want to process between progress reports?
    int blen; // Length of sequence to trim off from start to finish.
    int offset; // Number of bases at the start of the inline barcodes to skip for low quality.
    char ***rescaler; // Three-dimensional rescaler array. Size: [readlen, 39, 4] (length of reads, number of original quality scores, number of bases)
    char *rescaler_path; // Path to rescaler for
} mssi_settings_t;

typedef struct mark_splitter {
    FILE **tmp_out_handles_r1;
    FILE **tmp_out_handles_r2;
    int n_nucs;
    int n_handles;
    char **fnames_r1;
    char **fnames_r2;
} mark_splitter_t;


typedef struct sort_overlord {
    mark_splitter_t *splitter;
    //FILE **sort_out_handles_r1;
    //FILE **sort_out_handles_r2;
    char **out_fnames_r1;
    char **out_fnames_r2;
} sort_overlord_t;


char *barcode_mem_view(kseq_t *seq);
char ***parse_rescaler(char *qual_rescale_fname);
void free_mssi_settings(mssi_settings_t settings);

//Multiply a phred score by this to convert a -10log_10(x) to a -2log_e(x)
#define LOG10E_X5_INV 0.460517018598809136803598290936872841520220297725754595206665580193514521935470496
//such as in the following macro:
#define LOG10_TO_CHI2(x) (x) * LOG10E_X5_INV

// Throws an error when inferring barcode length.
#define MAX_BARCODE_LENGTH 30


// Calls incomplete gamma complement from CEPHES.
extern double igamc(double x, double y);


// Converts a chi2 sum into a p value.
static inline double igamc_pvalues(int num_pvalues, double x)
{
    if(x < 0) {
        return 1.0;
    }
    else {
        return igamc((double)num_pvalues, x / 2.0);
    }
}


//Converts a nucleotide in a char * into an index for the phred_sums and nuc_counts arrays.
#define NUC_TO_POS(character, nuc_indices)                                   \
    switch(character) {                                                      \
        case 'A': nuc_indices[0] = 0; nuc_indices[1] = 0; break;             \
        case 'C': nuc_indices[0] = 1; nuc_indices[1] = 1; break;             \
        case 'G': nuc_indices[0] = 2; nuc_indices[1] = 2; break;             \
        case 'T': nuc_indices[0] = 3; nuc_indices[1] = 3; break;             \
        default: nuc_indices[0] = 0; nuc_indices[1] = 4; break;              \
    }


typedef struct KingFisher {
    int **nuc_counts; // Count of nucleotides of this form
    double **phred_sums; // Sums of -10log10(p-value)
    int length; // Number of reads in family
    int readlen; // Length of reads
    char *max_phreds; // Maximum phred score observed at position. Use this as the final sequence for the quality to maintain compatibility with GATK and other tools.
    char *barcode;
    char pass_fail;
} KingFisher_t;



inline KingFisher_t init_kf(int readlen) {
    int **nuc_counts = (int **)malloc(readlen * sizeof(int *));
    double **phred_sums = (double **)malloc(sizeof(double *) * readlen);
    for(int i = 0; i < readlen; i++) {
        nuc_counts[i] = (int *)calloc(5, sizeof(int)); // One each for A, C, G, T, and N
        phred_sums[i] = (double *)calloc(4, sizeof(double)); // One for each nucleotide
    }
    KingFisher_t fisher = {
        .nuc_counts = nuc_counts,
        .phred_sums = phred_sums,
        .length = 0,
        .readlen = readlen,
        .max_phreds = (char *)calloc(readlen + 1, 1) // Keep track of the maximum phred score observed at position.
    };
    return fisher;
}


inline void destroy_kf(KingFisher_t *kfp) {
    fprintf(stderr, "Starting to destroy kfp with readlen %i.\n", kfp->readlen);
    for(int i = 0; i < kfp->readlen; ++i) {
        fprintf(stderr, "Starting to destroy.\n");
        fprintf(stderr, "Freeing nuc_counts and phred_sums %i.", i);
        free(kfp->nuc_counts[i]);
        free(kfp->phred_sums[i]);
    }
    free(kfp->nuc_counts);
    free(kfp->phred_sums);
    free(kfp->max_phreds);
    free(kfp->barcode);
}


inline void clear_kf(KingFisher_t *kfp) {
    for(int i = 0; i < kfp->readlen; i++) {
        memset(kfp->nuc_counts[i], 0, 5 * sizeof(int)); // And these.
        memset(kfp->phred_sums[i], 0, 4 * sizeof(double)); // Sets these to 0.
    }
    memset(kfp->max_phreds, 0, kfp->readlen); //Turn it back into an array of nulls.
    kfp->length = 0;
    return;
}

inline void pushback_kseq(KingFisher_t *kfp, kseq_t *seq, int *nuc_indices, int blen) {
    fprintf(stderr, "Pushing back kseq with read length %i\n", kfp->readlen);
    for(int i = 0; i < kfp->readlen; i++) {
        NUC_TO_POS((seq->seq.s[i]), nuc_indices);
        kfp->nuc_counts[i][nuc_indices[0]] += 1;
        kfp->phred_sums[i][nuc_indices[1]] += seq->qual.s[i] - 33;
        if(seq->qual.s[i] > kfp->max_phreds[i]) {
            kfp->max_phreds[i] = seq->qual.s[i];
        }
    }
    if(kfp->length == 0) {
        char *bs_ptr = barcode_mem_view(seq);
        kfp->pass_fail = (char)*(bs_ptr- 5);
        kfp->barcode = (char *)calloc(blen + 1, sizeof(char));
        memcpy(kfp->barcode, bs_ptr, blen);
    }
    kfp->length++; // Increment
    fprintf(stderr, "New length of kfp: %i. BTW, readlen for kfp: %i.\n", kfp->length, kfp->readlen);
    return;
}

/*
 * Warning: returns a NULL upon not finding a second pipe symbol.
 * This is *NOT* a properly null-terminated string.
 */
char *barcode_mem_view(kseq_t *seq) {
    int hits = 0;
    for(int i = 0; i < seq->comment.l; i++) {
        if(seq->comment.s[i] == '|') {
            if(!hits) {
                hits += 1;
            }
            else {
                return (char *)(seq->comment.s + i + 4); // 4 for "|BS="
            }
        }
    }
    return NULL;
}


inline int ARRG_MAX(KingFisher_t *kfp, int index) {
    /*
    fprintf(stderr, "Current values of phred_sums: %f,%f,%f,%f\n",
                    (double)kfp->phred_sums[index][0],
                    (double)kfp->phred_sums[index][1],
                    (double)kfp->phred_sums[index][2],
                    (double)kfp->phred_sums[index][3]);
    */
    if(kfp->phred_sums[index][3] > kfp->phred_sums[index][2] &&
       kfp->phred_sums[index][3] > kfp->phred_sums[index][1] &&
       kfp->phred_sums[index][3] > kfp->phred_sums[index][0]) {
        return 3;
    }
    else if(kfp->phred_sums[index][2] > kfp->phred_sums[index][1] &&
            kfp->phred_sums[index][2] > kfp->phred_sums[index][0]) {
        return 2;
    }
    else if(kfp->phred_sums[index][1] > kfp->phred_sums[index][0]) {
        return 1;
    }
    else {
        return 0;
    }
}

char ARRG_MAX_TO_NUC(int argmaxret) {
    switch (argmaxret) {
        case 1: return 'C';
        case 2: return 'G';
        case 3: return 'T';
        default: return 'A';
    }
}

inline int pvalue_to_phred(double pvalue) {
    return (int)(-10 * log10(pvalue));
}

inline void fill_csv_buffer(int readlen, int *arr, char *buffer, char *prefix) {
    char tmpbuf[10];
    sprintf(buffer, prefix);
    for(int i = 0; i < readlen; i++) {
        sprintf(tmpbuf, ",%i", arr[i]);
        strcat(buffer, tmpbuf);
    }
}

inline void fill_pv_buffer(KingFisher_t *kfp, int *phred_values, char *buffer) {
    fill_csv_buffer(kfp->readlen, phred_values, buffer, "PV:B:");
    return;
}

inline void fill_fa_buffer(KingFisher_t *kfp, int *agrees, char *buffer) {
    fill_csv_buffer(kfp->readlen, agrees, buffer, "FA:B:");
    return;
}

static inline void dmp_process_write(KingFisher_t *kfp, FILE *handle, int blen) {
    char pass_fail;
    char name_buffer[120];
    //1. Argmax on the phred_sums arrays, using that to fill in the new seq and
    char *cons_seq = (char *)malloc((kfp->readlen + 1) * sizeof(char));
    //buffer[0] = '@'; Set this later?
    int *cons_quals = (int *)malloc((kfp->readlen) * sizeof(int));
    int *agrees = (int *)malloc((kfp->readlen) * sizeof(int));
    cons_seq[kfp->readlen] = '\0'; // Null-terminate it.
    int argmaxret;
    for(int i = 0; i < kfp->readlen; i++) {
        argmaxret = ARRG_MAX(kfp, i);
        cons_seq[i] = ARRG_MAX_TO_NUC(argmaxret);
        cons_quals[i] = pvalue_to_phred(igamc_pvalues(kfp->length, LOG10_TO_CHI2((kfp->phred_sums[i][argmaxret]))));
        agrees[i] = kfp->nuc_counts[i][argmaxret];
    }
    cons_seq[kfp->readlen] = '\0'; // Null-terminal cons_seq.
    char FABuffer[1000];
    fill_fa_buffer(kfp, agrees, FABuffer);
    //fprintf(stderr, "FA buffer: %s.\n", FABuffer);
    char PVBuffer[1000];
    fill_pv_buffer(kfp, cons_quals, PVBuffer);
    char FPBuffer[7];
    sprintf(FPBuffer, "FP:i:%c", kfp->pass_fail);
    name_buffer[0] = '@';
    memcpy((char *)(name_buffer + 1), kfp->barcode, blen);
    name_buffer[1 + blen] = '\0';
    //fprintf(stderr, "Name buffer: %s\n", name_buffer);
    char arr_tag_buffer[2000];
    sprintf(arr_tag_buffer, "%s\t%s\t%s\n%s\n+\n%s\n", FABuffer, PVBuffer, FPBuffer, cons_seq, kfp->max_phreds);
    //fprintf(stderr, "Output result: %s %s", name_buffer, arr_tag_buffer);
    fprintf(handle, "%s %s", name_buffer, arr_tag_buffer);
    return;
}


#ifndef METASYNTACTIC_FNAME_BUFLEN
#define METASYNTACTIC_FNAME_BUFLEN 100
#endif

int ipow(int base, int exp);
void FREE_SPLITTER(mark_splitter_t var);


#ifndef KSEQ_2_FQ_INLINE
#define KSEQ_2_FQ_INLINE(handle, read, barcode, pass_fail, tmp_n_str, readlen, n_len) memcpy(tmp_n_str, read->seq.s, n_len * sizeof(char));\
        memset(tmp_n_str, 78, n_len);\
        fprintf(handle, "@%s ~#!#~|FP=%c|BS=%s\n%s\n+\n%s\n",\
                read->name.s, pass_fail, barcode, tmp_n_str, read->qual.s)
#endif


#ifndef KSEQ_2_FQ
#define KSEQ_2_FQ(handle, read, index, pass_fail) fprintf(handle, \
        "@%s ~#!#~|FP=%c|BS=%s\n%s\n+\n%s\n",\
    read->name.s, pass_fail, index->seq.s, read->seq.s, read->qual.s)
#endif


#ifndef FREE_SETTINGS
#define FREE_SETTINGS(settings) free(settings.output_basename);\
    free(settings.index_fq_path);\
    free(settings.input_r1_path);\
    free(settings.input_r2_path);
#endif



inline void free_mssi_settings(mssi_settings_t settings) {
    free(settings.output_basename);
    free(settings.input_r1_path);
    free(settings.input_r2_path);
    if(settings.rescaler_path) free(settings.rescaler_path);
    if(settings.rescaler) free(settings.rescaler);
    return;
}

#ifndef FREE_MSSI_SETTINGS
#define FREE_MSSI_SETTINGS(settings) free(settings.output_basename);\
    free(settings.input_r1_path);\
    free(settings.input_r2_path);\
    if(settings.rescaler_path) free(settings.rescaler_path);\
    if(settings.rescaler) free(settings.rescaler);
#endif


inline char test_hp_inline(char *barcode, int length, int threshold)
{
    int run = 0;
    char last = '\0';
    for(int i = 0; i < length; i++){
        if(barcode[i] == 'N') {
            return '0';
        }
        if(barcode[i] == last) {
            run += 1;
        }
        else {
            run = 0;
            last = barcode[i];
        }
    }
    return (run < threshold) ? '1': '0';
}

inline char test_hp(kseq_t *seq, int threshold)
{
    int run = 0;
    char last = '\0';
    for(int i = 0; i < seq->seq.l; i++){
        if(seq->seq.s[i] == 'N') {
            return '0';
        }
        if(seq->seq.s[i] == last) {
            run += 1;
        }
        else {
            run = 0;
            last = seq->seq.s[i];
        }
    }
    return (run < threshold) ? '1': '0';
}


// Functions
inline int ipow(int base, int exp)
{
    int result = 1;
    while (exp)
    {
        if (exp & 1)
            result *= base;
        exp >>= 1;
        base *= base;
    }

    return result;
}


void FREE_SPLITTER(mark_splitter_t var){
    for(int i = 0; i < var.n_handles; i++) {
        fclose(var.tmp_out_handles_r1[i]);
        fclose(var.tmp_out_handles_r2[i]);
        free(var.fnames_r1[i]);
        free(var.fnames_r2[i]);
    }
    free(var.tmp_out_handles_r1);
    free(var.tmp_out_handles_r2);
    return;
}

#define char_to_num(character, increment) switch(character) {\
    case 'C' : increment = 1; break;\
    case 'G' : increment = 2; break;\
    case 'T' : increment = 3; break;\
    default: increment = 0; break;\
    }


inline uint64_t get_binnerl(char *barcode, int length) {
    uint64_t bin = 0;
    int inc_binner;
    size_t count = 0;
    for(int i = length;i > 0;i--){
        char_to_num(barcode[i - 1], inc_binner);
        bin += (ipow(4, count) * inc_binner);
        count++;
    }
    return bin;
}


inline int get_binner(char *barcode, int length) {
    int bin = 0;
    int inc_binner;
    size_t count = 0;
    for(int i = length;i > 0;i--){
        char_to_num(barcode[i - 1], inc_binner);
        bin += (ipow(4, count) * inc_binner);
        count++;
    }
    return bin;
}


mark_splitter_t init_splitter(mss_settings_t* settings_ptr)
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


static void splitmark_core(kseq_t *seq1, kseq_t *seq2, kseq_t *seq_index,
                    mss_settings_t settings, mark_splitter_t splitter)
{
    int l1, l2, l_index, bin;
    int count = 0;
    char pass_fail;
    while ((l1 = kseq_read(seq1)) >= 0) {
        count += 1;
        if(!(count % settings.notification_interval)) {
            fprintf(stderr, "Number of records processed: %i.\n", count);
        }
        // Iterate through second fastq file.
        l_index = kseq_read(seq_index);
        l2 = kseq_read(seq2);
        if (l_index < 0) {
            fprintf(stderr, "Index return value for kseq_read less than "
                            "0. Are the fastqs different sizes? Abort!\n");
            FREE_SETTINGS(settings);
            FREE_SPLITTER(splitter);
            exit(EXIT_FAILURE);
        }
        if (l2 < 0) {
            fprintf(stderr, "Read 2 return value for kseq_read less than "
                            "0. Are the fastqs different sizes? Abort!\n");
            FREE_SETTINGS(settings);
            FREE_SPLITTER(splitter);
            exit(EXIT_FAILURE);
        }
        pass_fail = test_hp(seq_index, settings.hp_threshold);
        //fprintf(stdout, "Randomly testing to see if the reading is working. %s", seq1->seq.s);
        bin = get_binner(seq_index->seq.s, settings.n_nucs);
        KSEQ_2_FQ(splitter.tmp_out_handles_r1[bin], seq1, seq_index, pass_fail);
        KSEQ_2_FQ(splitter.tmp_out_handles_r2[bin], seq2, seq_index, pass_fail);
    }
}


inline int infer_barcode_length(char *bs_ptr) {
    int ret = 0;
    for (;;ret++) {
        if(bs_ptr[ret] == '\0') return ret;
#if !NDEBUG
        if(ret > MAX_BARCODE_LENGTH) {
            fprintf(stderr, "Inferred barcode length greater than max (%i). Abort!\n", MAX_BARCODE_LENGTH);
            exit(1);
        }
#endif
    }
}

#define rescaling_test(settings_ptr) (settings_ptr->rescaler_path)

inline int rescale_qscore(int qscore, int cycle, char base, char ***rescaler) {
    return rescaler[cycle][qscore - 2][base];
}


inline void pushback_rescaled_kseq(KingFisher_t *kfp, kseq_t *seq, char ***rescaler, int *nuc_indices, int blen) {
    fprintf(stderr, "Pushing back kseq with read length %i\n", kfp->readlen);
    for(int i = 0; i < kfp->readlen; i++) {
        NUC_TO_POS((seq->seq.s[i]), nuc_indices);
        kfp->nuc_counts[i][nuc_indices[0]] += 1;
        kfp->phred_sums[i][nuc_indices[1]] += rescale_qscore(seq->qual.s[i] - 33, i, seq->seq.s[i], rescaler);
        if(seq->qual.s[i] > kfp->max_phreds[i]) {
            kfp->max_phreds[i] = seq->qual.s[i];
        }
    }
    if(kfp->length == 0) {
        char *bs_ptr = barcode_mem_view(seq);
        kfp->pass_fail = (char)*(bs_ptr- 5);
        kfp->barcode = (char *)calloc(blen + 1, sizeof(char));
        memcpy(kfp->barcode, bs_ptr, blen);
    }
    kfp->length++; // Increment
    fprintf(stderr, "New length of kfp: %i. BTW, readlen for kfp: %i.\n", kfp->length, kfp->readlen);
    return;
}

inline char ***parse_rescaler(char *qual_rescale_fname) {
    char ***omgz = NULL;
    fprintf(stderr, "raise NotImplementedError('Hey, where do you think you're going? You're going NOWHERE. I gotchu foh 3 minutes. 3 minutes of playtime!')\n");
    exit(137);
    return omgz;
}


// mseq is a mutable struct holding kseq's information.

typedef struct mseq {
    char *name;
    char *comment;
    char *seq;
    char *qual;
    char *barcode;
    int l;
    int blen;
} mseq_t;


#define NUC_CMPL(character, ret) \
    switch(character) {\
        case 'A': ret = 'T'; break;\
        case 'C': ret = 'G'; break;\
        case 'G': ret = 'C'; break;\
        case 'T': ret = 'A'; break;\
        default: ret = 'N'; break;\
    }

inline int nuc_cmp(char forward, char reverse) {
    char tmpchar;
    NUC_CMPL(reverse, tmpchar)
    return forward - reverse;
}


inline int crc_flip(mseq_t *mvar, char *barcode, int blen, int readlen) {
    int cmp_ret;
    for(int i = 0; i < blen; ++i) {
        cmp_ret = nuc_cmp(barcode[i], barcode[blen - i - 1]);
        if(cmp_ret < 0) {
            return 0; // It's lexicographically lower as is. Don't flip!
        }
        else if(cmp_ret > 0) {
            return 1; // It's not lexicographically lower as it is. Flip!
        }
    }
    for(int i = 0; i < readlen; ++i) {
        cmp_ret = nuc_cmp(barcode[i], barcode[readlen - i - 1]);
        if(cmp_ret < 0) {
            return 0; // It's lexicographically lower as is. Don't flip!
        }
        else if(cmp_ret > 0) {
            return 1; // It's not lexicographically lower as it is. Flip!
        }
    }
    return 0; // Both barcode and read are lexicographically identical forward and reverse... is that possible? Eh. Don't flip.
}

typedef struct tmp_mseq {
    char *tmp_seq;
    char *tmp_qual;
    char *tmp_barcode;
    int readlen;
    int blen;
} tmp_mseq_t;

tmp_mseq_t init_tmp_mseq(int readlen, int blen) {
	char *tmp_seq = (char *)malloc(readlen);
	char *tmp_qual = (char *)malloc(readlen);
	char *tmp_barcode = (char *)malloc(blen);
    tmp_mseq_t ret = {
        .tmp_seq = tmp_seq,
        .tmp_qual = tmp_qual,
        .tmp_barcode = tmp_barcode,
        .readlen = readlen,
        .blen = blen
    };
    return ret;
}


void destroy_tmp_mseq(tmp_mseq_t mvar) {
    free(mvar.tmp_seq);
    free(mvar.tmp_qual);
    free(mvar.tmp_barcode);
    mvar.readlen = 0;
    mvar.blen = 0;
}


inline void crc_mseq(mseq_t *mvar, tmp_mseq_t *tmp)
{
    if(!crc_flip(mvar, mvar->barcode, tmp->blen, tmp->readlen)) return;
    for(int i = 0; i < tmp->readlen; i++) {
        NUC_CMPL(mvar->seq[tmp->readlen - i - 1], tmp->tmp_seq[i]);
        tmp->tmp_qual[i] = mvar->qual[tmp->readlen - i - 1];
    }
    for(int i = 0; i < tmp->blen; i++) {
        NUC_CMPL(mvar->seq[tmp->blen - i - 1], tmp->tmp_barcode[i]);
    }
    memcpy(mvar->qual, tmp->tmp_qual, tmp->readlen * sizeof(char));
    memcpy(mvar->seq, tmp->tmp_seq, tmp->readlen * sizeof(char));
    memcpy(mvar->barcode, tmp->tmp_barcode, tmp->blen * sizeof(char));
    return;
}


inline void mseq_rescale_init(kseq_t *seq, mseq_t *ret, char ***rescaler, tmp_mseq_t *tmp) {
    if(!seq) {
        ret = NULL;
        return;
    }
    ret->name = strdup(seq->name.s);
    ret->comment = strdup(seq->comment.s);
    ret->seq = strdup(seq->seq.s);
    if(!rescaler) ret->qual = strdup(seq->qual.s);
    else {
        ret->qual = (char *)malloc((seq->seq.l + 1) * sizeof(char));
        ret->qual[seq->seq.l] = '\0'; // Null-terminate this string.
        for(int i = 0; i < seq->seq.l; i++) {
            ret->qual[i] = rescale_qscore(ret->qual[i], i, ret->seq[i], rescaler);
        }
    }
    ret->l = seq->seq.l;
    crc_mseq(ret, tmp);
}


inline mseq_t init_rescale_revcmp_mseq(kseq_t *seq, char *barcode, char ***rescaler, tmp_mseq_t *tmp) {
    mseq_t ret = {
            .name = NULL,
            .comment = NULL,
            .seq = NULL,
            .qual = NULL,
			.barcode = barcode, // barcode still belongs to the argument variable!
            .l = 0,
			.blen = 0
    };
    mseq_rescale_init(seq, &ret, rescaler, tmp);
    return ret;
}


inline void mseq_free(mseq_t *mvar) {
    free(mvar->name);
    free(mvar->comment);
    free(mvar->seq);
    free(mvar->qual);
    mvar->l = 0;
    mvar->blen = 0;
    return;
}
