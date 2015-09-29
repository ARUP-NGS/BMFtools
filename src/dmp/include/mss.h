#pragma once

#include "stdio.h"
#include "mssi.h"
#include "binner.h"


#ifndef MAX_BARCODE_LENGTH
#define MAX_BARCODE_LENGTH 30
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
    free(settings.input_r2_path)
#endif

#ifndef METASYNTACTIC_FNAME_BUFLEN
#define METASYNTACTIC_FNAME_BUFLEN 100
#endif


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


typedef struct mark_splitter {
    FILE **tmp_out_handles_r1;
    FILE **tmp_out_handles_r2;
    int n_nucs;
    int n_handles;
    char **fnames_r1;
    char **fnames_r2;
} mark_splitter_t;

mark_splitter_t init_splitter(mss_settings_t* settings_ptr);

static inline void FREE_SPLITTER(mark_splitter_t var)
{
    for(int i = 0; i < var.n_handles; i++) {
        //fprintf(stderr, "Now trying to close file #%i with filename %s.\n", i, var.fnames_r1[i]);
        //fprintf(stderr, "Now trying to access FILE * with number %i.\n", i);
        if(var.fnames_r1[i]) {
            free(var.fnames_r1[i]);
            var.fnames_r1[i] = NULL;
        }
        if(var.fnames_r2[i]) {
            free(var.fnames_r2[i]);
            var.fnames_r2[i] = NULL;
        }
    }
    free(var.tmp_out_handles_r1);
    free(var.tmp_out_handles_r2);
    return;
}


static void splitmark_core(kseq_t *seq1, kseq_t *seq2, kseq_t *seq_index,
                           mss_settings_t settings, mark_splitter_t splitter)
{
    int l1, l2, l_index, bin;
    int count = 1;
    char pass_fail;
    while ((l1 = kseq_read(seq1)) >= 0 && (l2 = kseq_read(seq2) >= 0)
            && (l_index = kseq_read(seq_index)) >= 0) {
        if(!(++count % settings.notification_interval)) {
            fprintf(stderr, "Number of records processed: %i.\n", count);
        }
        // Iterate through second fastq file.
        pass_fail = test_hp(seq_index, settings.hp_threshold);
        //fprintf(stdout, "Randomly testing to see if the reading is working. %s", seq1->seq.s);
        bin = get_binner(seq_index->seq.s, settings.n_nucs);
        KSEQ_2_FQ(splitter.tmp_out_handles_r1[bin], seq1, seq_index, pass_fail);
        KSEQ_2_FQ(splitter.tmp_out_handles_r2[bin], seq2, seq_index, pass_fail);
    }
}

static void splitmark_core1(char *r1fq, char *r2fq, char *index_fq,
                            mss_settings_t *settings)
{
    gzFile fp_read1, fp_read2, fp_index;
    kseq_t *seq1, *seq2, *seq_index;
    mark_splitter_t splitter = init_splitter(settings);
    fp_read1 = gzopen(r1fq, "r");
    fp_read2 = gzopen(r2fq, "r");
    fp_index = gzopen(index_fq, "r");
    seq1 = kseq_init(fp_read1);
    seq2 = kseq_init(fp_read2);
    seq_index = kseq_init(fp_index);
    int l1, l2, l_index, bin;
    int count = 0;
    char pass_fail;
    while ((l1 = kseq_read(seq1)) >= 0 && (l2 = kseq_read(seq2) >= 0)
            && (l_index = kseq_read(seq_index)) >= 0) {
        if(!(count++ % settings->notification_interval)) {
            fprintf(stderr, "Number of records processed: %i.\n", count);
        }
        // Iterate through second fastq file.
        pass_fail = test_hp(seq_index, settings->hp_threshold);
        //fprintf(stdout, "Randomly testing to see if the reading is working. %s", seq1->seq.s);
        bin = get_binner(seq_index->seq.s, settings->n_nucs);
        KSEQ_2_FQ(splitter.tmp_out_handles_r1[bin], seq1, seq_index, pass_fail);
        KSEQ_2_FQ(splitter.tmp_out_handles_r2[bin], seq2, seq_index, pass_fail);
    }
    kseq_destroy(seq1);
    kseq_destroy(seq2);
    kseq_destroy(seq_index);
    gzclose(fp_read1);
    gzclose(fp_read2);
    gzclose(fp_index);
    return;
}


inline int infer_barcode_length(char *bs_ptr)
{
    int ret = 0;
    for (;;++ret) {
        if(bs_ptr[ret] == '\0' || bs_ptr[ret] == '|') {
#if !NDEBUG
            fprintf(stderr, "Returning blen as %i.\n", ret);
#endif
            return ret;
        }
#if !NDEBUG
        fprintf(stderr, "Length: %i. Char: %c.\n", ret, bs_ptr[ret]);
        if(ret > MAX_BARCODE_LENGTH) {
            fprintf(stderr, "Inferred barcode length greater than max (%i). Abort!\n", MAX_BARCODE_LENGTH);
            exit(1);
        }
#endif
    }
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
