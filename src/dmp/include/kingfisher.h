#pragma once

#include "stdio.h"
#include "math.h"
#include "charcmp.h"
#include "khash.h"
#include "uthash.h"
#include "khash.h"

#ifndef MAX_BARCODE_LENGTH
#define MAX_BARCODE_LENGTH 30
#endif

int nuc2num(char character);


// Memory costs
// max_phreds --> 1 * readlen + 8 (ptr)
// barcode --> [MAX_BARCODE_LENGTH + 1] + 8 (ptr)
// pass_fail --> 1
// nuc_counts --> readlen * 8 (ptrs) + 8 (ptr) + (readlen * 5 * sizeof(int)) [20] --> 6 * readlen + 1
// (nuc_counts = 28 * readlen + 8
// phred_sums --> readlen * 8 (ptrs) + 8 (ptr) + (readlen * 4 * sizeof(double)) [40 * readlen + 8]
// readlen --> 4
// length --> 1
// 69 * readlen + 46
// 4900 + 49 --> 5kB per barcode

typedef struct KingFisher {
    int **nuc_counts; // Count of nucleotides of this form
    double **phred_sums; // Sums of -10log10(p-value)
    int length; // Number of reads in family
    int readlen; // Length of reads
    char *max_phreds; // Maximum phred score observed at position. Use this as the final sequence for the quality to maintain compatibility with GATK and other tools.
    char barcode[MAX_BARCODE_LENGTH + 1];
    char pass_fail;
    int n_rc;
} KingFisher_t;



typedef struct tmpbuffers {
    char name_buffer[120];
    char PVBuffer[1000];
    char FABuffer[1000];
    char cons_seq_buffer[300];
    int cons_quals[300];
    int agrees[300];
} tmpbuffers_t;

extern double igamc(double x, double y);


//Multiply a phred score by this to convert a -10log_10(x) to a -2log_e(x)
#define LOG10E_X5_INV 0.460517018598809136803598290936872841520220297725754595206665580193514521935470496
//such as in the following macro:
#define LOG10_TO_CHI2(x) (x) * LOG10E_X5_INV



inline int pvalue_to_phred(double pvalue)
{
    return (int)(-10 * log10(pvalue));
}

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


inline KingFisher_t init_kf(int readlen)
{
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
        .max_phreds = (char *)calloc(readlen + 1, 1), // Keep track of the maximum phred score observed at position.
        .n_rc = 0,
        .pass_fail = '1'
    };
    return fisher;
}


inline void destroy_kf(KingFisher_t *kfp)
{
#if !NDEBUG
    fprintf(stderr, "Starting to destroy kfp with readlen %i.\n", kfp->readlen);
#endif
    for(int i = 0; i < kfp->readlen; ++i) {
#if !NDEBUG
        fprintf(stderr, "Starting to destroy.\n");
        fprintf(stderr, "Freeing nuc_counts and phred_sums %i.", i);
#endif
        free(kfp->nuc_counts[i]);
        free(kfp->phred_sums[i]);
    }
    free(kfp->nuc_counts);
    free(kfp->phred_sums);
    free(kfp->max_phreds);
}


inline void clear_kf(KingFisher_t *kfp)
{
    for(int i = 0; i < kfp->readlen; i++) {
        memset(kfp->nuc_counts[i], 0, 5 * sizeof(int)); // And these.
        memset(kfp->phred_sums[i], 0, 4 * sizeof(double)); // Sets these to 0.
    }
    memset(kfp->max_phreds, 0, kfp->readlen); //Turn it back into an array of nulls.
    kfp->length = 0;
    return;
}


inline int ARRG_MAX(KingFisher_t *kfp, int index)
{
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

inline char ARRG_MAX_TO_NUC(int argmaxret)
{
    switch (argmaxret) {
        case 1: return 'C';
        case 2: return 'G';
        case 3: return 'T';
        default: return 'A';
    }
}


inline void fill_csv_buffer(int readlen, int *arr, char *buffer, char *prefix, char typecode)
{
    char tmpbuf[20];
    sprintf(buffer, "%s%c", prefix, typecode);
    for(int i = 0; i < readlen; i++) {
        sprintf(tmpbuf, ",%i", arr[i]);
        strcat(buffer, tmpbuf);
    }
}

inline void fill_pv_buffer(KingFisher_t *kfp, int *phred_values, char *buffer)
{
    fill_csv_buffer(kfp->readlen, phred_values, buffer, "PV:B:", 'I');
    return;
}

inline void fill_fa_buffer(KingFisher_t *kfp, int *agrees, char *buffer)
{
    fill_csv_buffer(kfp->readlen, agrees, buffer, "FA:B:", 'I'); // Add in the "I" to type the array.
    return;
}


/*
 * TODO: Use tmpvals_t object to avoid allocating and deallocating each of these.
 */
static inline void dmp_process_write(KingFisher_t *kfp, FILE *handle, int blen, tmpbuffers_t *tmp)
{
#if !NDEBUG
    fprintf(stderr, "Hey, just got to dpw.\n");
#endif
    //1. Argmax on the phred_sums arrays, using that to fill in the new seq and
    //buffer[0] = '@'; Set this later?
    int argmaxret;
    tmp->cons_seq_buffer[kfp->readlen] = '\0'; // Null-terminal cons_seq.
    for(int i = 0; i < kfp->readlen; ++i) {
        argmaxret = ARRG_MAX(kfp, i);
        tmp->cons_seq_buffer[i] = ARRG_MAX_TO_NUC(argmaxret);
        tmp->cons_quals[i] = pvalue_to_phred(igamc_pvalues(kfp->length, LOG10_TO_CHI2((kfp->phred_sums[i][argmaxret]))));
        tmp->agrees[i] = kfp->nuc_counts[i][argmaxret];
    }
    fill_fa_buffer(kfp, tmp->agrees, tmp->FABuffer);
    //fprintf(stderr, "FA buffer: %s.\n", FABuffer);
    fill_pv_buffer(kfp, tmp->cons_quals, tmp->PVBuffer);
    tmp->name_buffer[0] = '@';
    memcpy((char *)(tmp->name_buffer + 1), kfp->barcode, blen);
    tmp->name_buffer[1 + blen] = '\0';
    //fprintf(stderr, "Name buffer: %s\n", tmp->name_buffer);
    //fprintf(stderr, "Output result: %s %s", tmp->name_buffer, arr_tag_buffer);
    fprintf(handle, "%s %s\t%s\tFP:i:%c\n%s\n+\n%s\n", tmp->name_buffer, tmp->FABuffer, tmp->PVBuffer, kfp->pass_fail, tmp->cons_seq_buffer, kfp->max_phreds);
    return;
}


inline int rescale_qscore(int readnum, int qscore, int cycle, char base, char ****rescaler)
{
    return rescaler[readnum][cycle][qscore - 2][nuc2num(base)];
}



void set_kf(int readlen, KingFisher_t ret)
{
    ret.length = 0;
    ret.readlen = readlen;
    ret.nuc_counts = (int **)malloc(readlen * sizeof(int *));
    ret.phred_sums = (double **)malloc(sizeof(double *) * readlen);
    ret.max_phreds = (char *)calloc(readlen + 1, sizeof(char)); // Keep track of the maximum phred score observed at position.
    for(int i = 0; i < readlen; i++) {
        ret.nuc_counts[i] = (int *)calloc(5, sizeof(int)); // One each for A, C, G, T, and N
        ret.phred_sums[i] = (double *)calloc(4, sizeof(double)); // One for each nucleotide
    }
    return;
}
