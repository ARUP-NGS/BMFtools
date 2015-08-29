#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <stddef.h>
#include <quadmath.h>
#include "include/kseq.h"

// Force declaration of all of kseq's types.
KSEQ_INIT(gzFile, gzread)

typedef __float128 float128_t;


#define LOG10E_X5_INV 0.4605170185988091368035982909368728415202202977257545952066655801935145219354704960471994410179196596683935568084572497266819050930165613513332L
//Multiply a phred score by this to convert a -10log_10(x) to a -2log_e(x)
//such as in the following macro:
#define LOG10_TO_CHI2(x) (x) * LOG10E_X5_INV

#define INV_CHI2_FROM_LOG10(log10int) -2 * log(1 - pow(10, log10int))
/*
 * Equivalent to the following, but type-general:
inline float128_t INV_CHI2_FROM_LOG10(int32_t log10int)
{
    return -2 * log(1 - pow(10, log10int));
}
*/
#define INV_LOG10_FROM_LOG10(log10int) -10 * log10(1 - pow(10, log10int))

#define MAX_BARCODE_LENGTH 36


extern float128_t igamcl(float128_t a, float128_t x);

// Converts a 
inline float128_t igamc_pvalues(int num_pvalues, float128_t x)
{
    if(x < 0) {
        return 1.0;
    }
    else {
#if !NDEBUG
        fprintf(stderr, "Now calling igamcl.\n");
#endif
        return igamcl(num_pvalues * 1., x / 2.0);
    }
}

#define NUC_TO_POS(character, nuc_indices)                            \
    switch(character) {                                               \
        case 'A': nuc_indices[0] = 0; nuc_indices[1] = 0; break;      \
        case 'C': nuc_indices[0] = 1; nuc_indices[1] = 1; break;      \
        case 'G': nuc_indices[0] = 2; nuc_indices[1] = 2; break;      \
        case 'T': nuc_indices[0] = 3; nuc_indices[1] = 3; break;      \
        default: nuc_indices[0] = 0; nuc_indices[1] = 4; break;       \
    }


/*
 * TODO: KingFisher finishing work.
 * A destructor for KingFisher.
 * FFF core.
 * Output result:
 * Rewrite the cFastFisherFlattening array work in C rather than Cython from MawCluster/BCFastq.pyx.
 * Use that array work to fill in the update_kf method.
 */

typedef struct KingFisher {
    int **nuc_counts; // Count of nucleotides of this form
    float128_t **phred_sums; // Sums of -10log10(p-value)
    int length; // Number of reads in family
    int readlen; // Length of reads
    char *max_phreds; // Maximum phred score observed at position. Use this as the final sequence for the quality to maintain compatibility with GATK and other tools.
} KingFisher_t;

inline KingFisher_t init_kf(size_t readlen) {
    int **nuc_counts = (int **)malloc(readlen * sizeof(int *));
    float128_t **phred_sums = (float128_t **)malloc(sizeof(float128_t *) * 4);
    for(int i = 0; i < readlen; i++) {
        nuc_counts[i] = (int *)calloc(5, sizeof(int)); // One each for A, C, G, T, and N
        phred_sums[i] = (float128_t *)calloc(4, sizeof(float128_t)); // One for each nucleotide
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
    for(int i = 0; i < kfp->readlen; i++) {
        free(kfp->nuc_counts[i]);
        free(kfp->phred_sums[i]);
    }
    free(kfp->nuc_counts);
    free(kfp->phred_sums);
    free(kfp->max_phreds);
}

inline void clear_kf(KingFisher_t *kfp) {
    for(int i = 0; i < kfp->readlen; i++) {
        memset(kfp->phred_sums[i], 0, 4 * sizeof(float128_t)); // Sets these to 0.
        memset(kfp->nuc_counts[i], 0, 5 * sizeof(int)); // And these.
    }
    memset(kfp->max_phreds, 0, kfp->readlen); //Turn it back into an array of nulls.
    kfp->length = 0;
    return;
}

/*
inline void update_nuc_counts(KingFisher_t *fisher, kseq_t *seq){

    fprintf(stderr, "update_kf for updating KingFisher_t is unimplemented. Abort!\n");
    exit(1);
    return;
}
*/

inline void pushback_kseq(KingFisher_t *kfp, kseq_t *seq, int *nuc_indices) {
    for(int i = 0; i < kfp->readlen; i++) {
        NUC_TO_POS((seq->seq.s[i]), nuc_indices);
        kfp->nuc_counts[i][nuc_indices[0]] += 1;
        for(int j = 0; j < 4; j++) {
            if(j == nuc_indices[1])
            kfp->phred_sums[i][j] += (seq->qual.s[i] - 33);
            else {
                kfp->phred_sums[i][j] += INV_LOG10_FROM_LOG10((seq->qual.s[i] - 33)); // Make sure we decrease our confidence in other base calls as well.
            }
        }
        if(seq->seq.s[i] > kfp->max_phreds[i]) {
            kfp->max_phreds[i] = seq->seq.s[i];
        }
    }
    kfp->length++; // Increment
    return;
}

/*
 * Warning: returns a NULL upon not finding a second pipe symbol.
 * This is *NOT* a properly null-terminated string.
 */
inline char *barcode_mem_view(kseq_t *seq) {
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

inline char ARRG_MAX_TO_NUC(int argmaxret) {
    switch (argmaxret) {
        case 1: return 'C';
        case 2: return 'G';
        case 3: return 'T';
        default: return 'A';
    }
}

inline int pvalue_to_phred(float128_t pvalue) {
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
    fill_csv_buffer(kfp->readlen, phred_values, buffer, "PV:i:");
    return;
}

inline void fill_fa_buffer(KingFisher_t *kfp, int *agrees, char *buffer) {
    fill_csv_buffer(kfp->readlen, agrees, buffer, "FM:i:");
    return;
}


inline void dmp_process_write(KingFisher_t *kfp, FILE *handle, char *bs_ptr, int blen) {
    int pass;
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
    char PVBuffer[1000];
    fill_pv_buffer(kfp, cons_quals, PVBuffer);
    pass = (int)*(bs_ptr - 5); // 5 for |BS=[ACTG]
#if !NDEBUG
    if(pass != 0 && pass != 1) {
        char buf[20];
        memcpy(buf, (bs_ptr - 5), 5 + blen + 1);
        buf[5 + blen + 1] = '\0';
        fprintf(stderr, buf);
        exit(1);
    }
#endif
    char FPBuffer[7];
    sprintf(FPBuffer, "FP:i:%i", pass);
    name_buffer[0] = '@';
    strncpy((char *)(name_buffer + 1), bs_ptr, blen);
    name_buffer[1 + blen] = '\0';
    char arr_tag_buffer[2000];
    sprintf(arr_tag_buffer, "%s\t%s\t%s\n%s\n+\n%s\n", FABuffer, PVBuffer, FPBuffer, cons_seq, kfp->max_phreds);
    fprintf(handle, "%s %s", name_buffer, arr_tag_buffer);
    //Make the strings to write to handle
    //Write the strings to handle
    return;
}
