#include <getopt.h>
#include <math.h>
#include <omp.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "include/mss.h" // contains structs and methods related to both inline and non-inline mss


typedef struct sort_overlord {
    mark_splitter_t *splitter;
    //FILE **sort_out_handles_r1;
    //FILE **sort_out_handles_r2;
    char **out_fnames_r1;
    char **out_fnames_r2;
} sort_overlord_t;


// Throws an error when inferring barcode length.
#define MAX_BARCODE_LENGTH 30


// Calls incomplete gamma complement from CEPHES.

char *barcode_mem_view(kseq_t *seq);
char ***parse_rescaler(char *qual_rescale_fname);


#ifndef METASYNTACTIC_FNAME_BUFLEN
#define METASYNTACTIC_FNAME_BUFLEN 100
#endif


int ipow(int base, int exp);


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


typedef struct outpost {
    khash_t(fisher) *hash;
    kseq_t *seq;
    char *bs_ptr;
    int blen;
    //int readlen; - Don't need readlen - already available as seq->seq.l;
    int *nuc_indices;
    khiter_t k;
} outpost_t;

/*
static inline void pushback_hash(outpost_t Navy)
{
    Navy.bs_ptr = barcode_mem_view(Navy.seq);
    int64_t bin = get_binnerl(Navy.bs_ptr, Navy.blen);
    fprintf(stderr, "Bin for pushing back hash: %i.", bin);
    Navy.k=kh_get(fisher, Navy.hash,
                  get_binner(Navy.bs_ptr, Navy.blen));
    if(Navy.k==kh_end(Navy.hash)) {
        fprintf(stderr, "New barcode! %s, %i.", Navy.bs_ptr, bin);
        kh_value(Navy.hash, Navy.k) = init_kf(Navy.seq->seq.l);
        pushback_kseq(&kh_value(Navy.hash, Navy.k), Navy.seq, Navy.nuc_indices, Navy.blen);
    }
    else {
        pushback_kseq(&kh_value(Navy.hash, Navy.k), Navy.seq, Navy.nuc_indices, Navy.blen);
    }
    return;
}
*/
