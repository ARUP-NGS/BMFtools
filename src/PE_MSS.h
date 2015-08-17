#include <zlib.h>
#include <stdio.h>
#include <stdint.h>
#include <regex.h>
#include <getopt.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include "include/kseq.h"
#include "include/sort/lh3sort.c"

#ifndef METASYNTACTIC_FNAME_BUFLEN
#define METASYNTACTIC_FNAME_BUFLEN 100
#endif
#ifndef RANDOM_STRING_LENGTH
#define RANDOM_STRING_LENGTH 30
#endif
#ifndef MAX_BARCODE_PREFIX_LENGTH
#define MAX_BARCODE_PREFIX_LENGTH 12
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
    char *tmp_split_basename;
} mss_settings_t;

typedef struct mark_splitter {
    FILE **tmp_out_handles_r1;
    FILE **tmp_out_handles_r2;
    int n_nucs;
    int n_handles;
    char **fnames_r1;
    char **fnames_r2;
} mark_splitter_t;


typedef struct sort_overlord {
    mark_splitter_t splitter;
    FILE **sort_out_handles_r1;
    FILE **sort_out_handles_r2;
    char **out_fnames_r1;
    char **out_fnames_r2;
} sort_overlord_t;


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

inline int lh3_sort_call(char *fname, char *outfname)
{
    int retvar;
    char **lh3_argv = (char **)malloc(6 * sizeof(char *));
    lh3_argv[1] = strdup("-t\'|\'");
    lh3_argv[2] = strdup("-k2,2");
    lh3_argv[3] = strdup("-o");
    lh3_argv[4] = strdup(outfname);
    lh3_argv[5] = strdup(fname);
    retvar = lh3_sort_main(6, lh3_argv);
    for(int i = 1; i < 6; i++) {
        free(lh3_argv[i]);
    }
    free(lh3_argv);
    return retvar;
}


inline void FREE_SPLITTER(mark_splitter_t var){
    for(int i = 0; i < var.n_handles; i++)
    {
        fclose(var.tmp_out_handles_r1[i]);
        fclose(var.tmp_out_handles_r2[i]);
        free(var.fnames_r1[i]);
        free(var.fnames_r2[i]);
    }
    free(var.tmp_out_handles_r1);
    free(var.tmp_out_handles_r2);
    return;
}

inline void apply_lh3_sorts(sort_overlord_t *dispatcher, mss_settings_t *settings)
{
    int abort = 0;
    int index = -1;
    omp_set_num_threads(settings->threads);
    #pragma omp parallel for
    for(int i = 0; i < dispatcher->splitter.n_handles; i++) {
        #pragma omp flush(abort)
        int ret = lh3_sort_call(dispatcher->splitter.fnames_r1[i], dispatcher->out_fnames_r1[i]);
        if(!ret) {
            abort = 1;
            index = i;
            #pragma omp flush (abort)
            #pragma omp flush (index)
        }
    }
    if(abort) {
        fprintf(stderr,
                "lh3 sort call failed for file handle %s. (Non-zero exit status). Abort!",
                dispatcher->splitter.fnames_r1[index]);
                //FREE_MP_SORTER(*dispatcher); // Delete allocated memory.
                // Will need to rewrite this for paired-end.
        exit(EXIT_FAILURE);
    }
    abort = 0;
    index = -1;
    #pragma omp parallel for
    for(int i = 0; i < dispatcher->splitter.n_handles; i++) {
        #pragma omp flush(abort)
        int ret = lh3_sort_call(dispatcher->splitter.fnames_r2[i], dispatcher->out_fnames_r2[i]);
        if(!ret) {
            abort = 1;
            index = i;
            #pragma omp flush (abort)
            #pragma omp flush (index)
        }
    }
    if(abort) {
        fprintf(stderr,
                "lh3 sort call failed for file handle %s. (Non-zero exit status). Abort!",
                dispatcher->splitter.fnames_r2[index]);
                //FREE_MP_SORTER(*dispatcher); // Delete allocated memory.
                // Will need to rewrite this for paired-end.
        exit(EXIT_FAILURE);
    }
    return;
}


#define char_to_num(character, increment) switch(character) {\
    case 'C' : increment = 1; break;\
    case 'G' : increment = 2; break;\
    case 'T' : increment = 3; break;\
    case 'N' : increment = -65535; break;\
    default: increment = 0; break;\
    }


inline int get_binner(char binner[], int length) {
    int bin = 0;
    int inc_binner;
    size_t count = 0;
    for(int i = length;i > 0;i--){
        char_to_num(binner[i - 1], inc_binner);
        bin += (ipow(4, count) * inc_binner);
        count++;
    }
    return bin;
}

inline mark_splitter_t init_splitter(mss_settings_t* settings_ptr)
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
        sprintf(tmp_buffer, "%s.tmp.%i.R1.fastq", settings_ptr->tmp_split_basename, i);
        ret.fnames_r1[i] = strdup(tmp_buffer);
        sprintf(tmp_buffer, "%s.tmp.%i.R2.fastq", settings_ptr->tmp_split_basename, i);
        ret.fnames_r2[i] = strdup(tmp_buffer);
        ret.tmp_out_handles_r1[i] = fopen(ret.fnames_r1[i], "w");
        ret.tmp_out_handles_r2[i] = fopen(ret.fnames_r2[i], "w");
    }
    return ret;
}
