#pragma once

#include "kseq_dec.h"
#include "cstr_utils.h"


#ifndef METASYNTACTIC_FNAME_BUFLEN
#define METASYNTACTIC_FNAME_BUFLEN 100
#endif
#ifndef MAX_N_BLENS
#define MAX_N_BLENS 10
#endif

typedef struct blens {
    int max_blen; // Last value in blens
    int min_blen; // Lowest value in blens
    int blens[MAX_N_BLENS]; // Array holding blens
    int n; // Number of blens to look for
    int current_blen;
} blens_t;


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
    int offset; // Number of bases at the start of the inline barcodes to skip for low quality.
    char *rescaler_path; // Path to rescaler for
    int threads;
    int run_hash_dmp;
    char *ffq_prefix; // Final fastq prefix
    char ****rescaler; // Four-dimensional rescaler array. Size: [readlen, 39, 4] (length of reads, number of original quality scores, number of bases)
    int blen;
    int blen1_2;
    int max_blen;
    int gzip_output;
} mssi_settings_t;

void free_mssi_settings(mssi_settings_t settings);


inline void free_mssi_settings(mssi_settings_t settings)
{
    free(settings.output_basename);
    free(settings.input_r1_path);
    free(settings.input_r2_path);
    if(settings.rescaler_path) free(settings.rescaler_path);
    return;
}

#ifndef FREE_MSSI_SETTINGS_PTR
#define FREE_MSSI_SETTINGS_PTR(settings) free(settings->output_basename);\
    free(settings->input_r1_path);\
    free(settings->input_r2_path);\
    if(settings->rescaler_path) free(settings->rescaler_path)
#endif

#ifndef FREE_MSSI_SETTINGS
#define FREE_MSSI_SETTINGS(settings) free(settings.output_basename);\
    free(settings.input_r1_path);\
    free(settings.input_r2_path);\
    if(settings.rescaler_path) free(settings.rescaler_path);
#endif


inline int test_homing_seq(kseq_t *seq1, kseq_t *seq2, mssi_settings_t *settings_ptr)
{
    if(!settings_ptr->homing_sequence) {
        return 1;
    }
    else {
        return memcmp(seq1->seq.s + (settings_ptr->blen / 2 + settings_ptr->offset),
                       settings_ptr->homing_sequence,
                       settings_ptr->homing_sequence_length) == 0;
    }
}


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


const char *crms_suffix = ".crms.split";

/*
 * Returns a null-terminated string with the default outfname.
 * Warning: Must be freed!
 */
static inline char *make_default_outfname(char *fname, const char *suffix)
{
    char buf[200];
    char *prefix = trim_ext(fname);
    strcpy(buf, prefix);
    strcat(buf, suffix);
    char *ret = strdup(buf);
    free(prefix);
    return ret;
}

static inline char *make_crms_outfname(char *fname)
{
    return make_default_outfname(fname, crms_suffix);
}
