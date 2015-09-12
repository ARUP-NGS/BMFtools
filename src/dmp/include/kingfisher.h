#pragma once

#include "stdio.h"
#include "math.h"
#include "charcmp.h"
#include "khash.h"


typedef struct KingFisher {
    int **nuc_counts; // Count of nucleotides of this form
    double **phred_sums; // Sums of -10log10(p-value)
    int length; // Number of reads in family
    int readlen; // Length of reads
    char *max_phreds; // Maximum phred score observed at position. Use this as the final sequence for the quality to maintain compatibility with GATK and other tools.
    char *barcode;
    char pass_fail;
} KingFisher_t;


KHASH_MAP_INIT_INT64(fisher, KingFisher_t) // Initialize a hashmap with int64 keys and KingFisher_t payload.
int nuc2num(char character);

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
        .max_phreds = (char *)calloc(readlen + 1, 1) // Keep track of the maximum phred score observed at position.
    };
    return fisher;
}


inline void destroy_kf(KingFisher_t *kfp)
{
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


inline void fill_csv_buffer(int readlen, int *arr, char *buffer, char *prefix)
{
    char tmpbuf[10];
    sprintf(buffer, prefix);
    for(int i = 0; i < readlen; i++) {
        sprintf(tmpbuf, ",%i", arr[i]);
        strcat(buffer, tmpbuf);
    }
}

inline void fill_pv_buffer(KingFisher_t *kfp, int *phred_values, char *buffer)
{
    fill_csv_buffer(kfp->readlen, phred_values, buffer, "PV:B:");
    return;
}

inline void fill_fa_buffer(KingFisher_t *kfp, int *agrees, char *buffer)
{
    fill_csv_buffer(kfp->readlen, agrees, buffer, "FA:B:");
    return;
}

static inline void dmp_process_write(KingFisher_t *kfp, FILE *handle, int blen)
{
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


inline int rescale_qscore(int qscore, int cycle, char base, char ***rescaler)
{
    return rescaler[cycle][qscore - 2][nuc2num(base)];
}
