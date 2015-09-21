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
#ifndef MAX_BARCODE_LENGTH
#define MAX_BARCODE_LENGTH 30
#endif

#ifndef kroundup32
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif


// Calls incomplete gamma complement from CEPHES.

char *barcode_mem_view(kseq_t *seq);


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


KHASH_MAP_INIT_INT64(fisher, KingFisher_t *) // Initialize a hashmap with uint64 keys and KingFisher_t payload.


static inline KingFisher_t *init_kfp(size_t readlen)
{
    KingFisher_t *ret = (KingFisher_t *)malloc(sizeof(KingFisher_t));
    ret->length = 0; // Check to see if this is necessary after calloc - I'm pretty sure not.
    ret->n_rc = 0;
    ret->readlen = readlen;
#if !NDEBUG
    fprintf(stderr, "New read length for new kfp: %i. Pointer: %p.", ret->readlen, ret);
#endif
    ret->max_phreds = (char *)malloc((readlen + 1) * sizeof(char)), // Keep track of the maximum phred score observed at position.
    ret->nuc_counts = (uint16_t **)malloc(readlen * sizeof(uint16_t *));
    ret->phred_sums = (uint16_t **)malloc(readlen * sizeof(uint16_t *));
    for(int i = 0; i < readlen; ++i) {
        ret->nuc_counts[i] = (uint16_t *)calloc(5, sizeof(uint16_t)); // One each for A, C, G, T, and N
        ret->phred_sums[i] = (uint16_t *)calloc(4, sizeof(uint16_t)); // One for each nucleotide
    }
    ret->pass_fail = '1';
    return ret;
}



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
