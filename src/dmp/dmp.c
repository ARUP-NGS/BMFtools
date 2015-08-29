#include <getopt.h>
#include <quadmath.h>
#include "dmp.h"

//Function definitions
float128_t igamc_pvalues(int num_pvalues, float128_t x);
KingFisher_t init_kf(size_t readlen);
void pushback_kseq(KingFisher_t *fisher, kseq_t *seq, int *nuc_indices);
int bmftools_dmp_core(kseq_t *seq, FILE *out_handle, int readlen, int blen);
int ARRG_MAX(KingFisher_t *kfp, int index);
char ARRG_MAX_TO_NUC(int argmaxret);
int pvalue_to_phred(float128_t pvalue);
void fill_fa_buffer(KingFisher_t *kfp, int *agrees, char *buffer);

void print_usage() {
    fprintf(stderr, "This isn't ready to do anything yet. Oops.\n");
}

void print_opt_err(char *argv[], char optarg[]) {
    fprintf(stderr, "Invalid argument %s. See usage.\n", optarg);
    print_usage();
    exit(1);
}

int main(int argc, char* argv[]) {
    int n = 0;
    double f = 1337.;
    int c;
    while ((c = getopt(argc, argv, "n:f:")) > -1) {
        switch(c) {
            case 'n': n = atoi(optarg); break;
            case 'f': f = atof(optarg); break;
            default: print_opt_err(argv, optarg);
        }
    }
#if !NDEBUG
    fprintf(stderr, "IGAMC_PVALUES: %.30Lf\n", (long double)igamc_pvalues((float128_t)n, (float128_t)f));
    fprintf(stderr, "degrees of freedom:  %i. Chi2 sum: %f. Chi2 sum as a string: %s.\n", n, f, argv[2]);
    float128_t *arr = (float128_t *)calloc(1, sizeof(float128_t));
    int wtf = arr[0];
    fprintf(stderr, "A zeroed 128-bit float has value %i\n", wtf);
#endif
    return 0;
}


int bmftools_dmp_wrapper(char *input_path, char *output_path,
                         int readlen, int blen) {
    /*
     * Set output_path to NULL to write to stdout.
     * Set input_path to "-" or NULL to read from stdin.
     */
    FILE *in_handle;
    FILE *out_handle;
    if(input_path[0] == '-' || !input_path) in_handle = stdin;
    else {
        in_handle = fopen(input_path, "r");
    }
    if(!output_path) out_handle = stdout;
    else {
        out_handle = fopen(output_path, "w");
    }
    gzFile fp = gzdopen(fileno(in_handle), "r");
    kseq_t *seq = kseq_init(fp);
    return bmftools_dmp_core(seq, out_handle, readlen, blen);
}


void dmp_process_write(KingFisher_t *kfp, FILE *handle, char *bs_ptr, int blen) {
	int pass;
	char name_buffer[120];
    //1. Argmax on the chi2sums arrays, using that to fill in the new seq and
    char *cons_seq = (char *)malloc((kfp->readlen + 1) * sizeof(char));
    //buffer[0] = '@'; Set this later?
    int *cons_quals = (int *)malloc((kfp->readlen) * sizeof(int));
    int *agrees = (int *)malloc((kfp->readlen) * sizeof(int));
    cons_seq[kfp->readlen] = '\0'; // Null-terminate it.
    int argmaxret;
    for(int i = 0; i < kfp->readlen; i++) {
        argmaxret = ARRG_MAX(kfp, i);
        cons_seq[i] = ARRG_MAX_TO_NUC(argmaxret);
        cons_quals[i] = pvalue_to_phred(igamc_pvalues(kfp->length, kfp->chi2sums[i][argmaxret]));
        agrees[i] = kfp->nuc_counts[i][argmaxret];
    }
    cons_seq[kfp->readlen] = '\0'; // Null-terminal cons_seq.
    char FABuffer[1000];
    fill_fa_buffer(kfp, agrees, FABuffer);
    char PVBuffer[1000];
    fill_pv_buffer(kfp, cons_quals, PVBuffer);
    pass = (int)*(bs_ptr - 2);
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


int bmftools_dmp_core(kseq_t *seq, FILE *out_handle, int readlen, int blen) {
    int l;
    int *nuc_indices = malloc(2 * sizeof(int));
    char *current_barcode = (char *)calloc(blen + 1, 1);
    KingFisher_t Holloway = init_kf(readlen);
    KingFisher_t *Hook = &Holloway;
    l = kseq_read(seq);
    if(l < 0) {
        fprintf(stderr, "Could not open input fastq. Abort!\n");
        exit(1);
    }
    char *bs_ptr = barcode_mem_view(seq); //Note: NOT null-terminated at the end of the barcode.
    memcpy(current_barcode, bs_ptr, blen);
    pushback_kseq(Hook, seq, nuc_indices); // Initialize Hook with
    while ((l = kseq_read(seq)) >= 0) {
        bs_ptr = barcode_mem_view(seq);
#if !NDEBUG
        if(!bs_ptr) { // bs_ptr is NULL
            fprintf(stderr, "Malformed fastq comment field - missing the second delimiter. Abort!\n");
            exit(1);
        }
#endif
        if(memcmp(bs_ptr, current_barcode, blen) == 0) {
            pushback_kseq(Hook, seq, nuc_indices);
        }
        else {
            dmp_process_write(Hook, out_handle, bs_ptr, blen);
            clear_kf(Hook); // Reset Holloway
            memcpy(current_barcode, bs_ptr, blen * sizeof(char)); // Update working barcode.
            pushback_kseq(Hook, seq, nuc_indices);
        }
    }
    if(Hook->length) { // If length is not 0
        dmp_process_write(Hook, out_handle, bs_ptr, blen);
    }
    destroy_kf(Hook);
    kseq_destroy(seq);

    return 0;
}
