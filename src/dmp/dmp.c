#include <getopt.h>
#include <quadmath.h>
#include "dmp.h"

//Function definitions
float128_t igamc_pvalues(int num_pvalues, float128_t x);
KingFisher_t init_kf(size_t readlen);
void pushback_kseq(KingFisher_t *fisher, kseq_t *seq, int *nuc_indices);
int bmftools_dmp_core(kseq_t *seq, FILE *out_handle, int blen);
int ARRG_MAX(KingFisher_t *kfp, int index);
char ARRG_MAX_TO_NUC(int argmaxret);
int pvalue_to_phred(float128_t pvalue);
void fill_fa_buffer(KingFisher_t *kfp, int *agrees, char *buffer);
void dmp_process_write(KingFisher_t *kfp, FILE *handle, char *bs_ptr, int blen);
int bmftools_dmp_wrapper(char *input_path, char *output_path,
                         int blen);

void print_usage(char *argv[]) {
    fprintf(stderr, "Usage: %s -l <blen> -o <output_path> <input_path>\n"
                    "Set output_path to \"-\" for stdout, "
                    "input_path to \"-\" for stdin.\n", argv[0]);
}

void print_opt_err(char *argv[], char optarg[]) {
    fprintf(stderr, "Invalid argument %s. See usage.\n", optarg);
    print_usage(argv);
    exit(1);
}

int main(int argc, char* argv[]) {
    int blen = 0;
#ifdef BMF_THREADS
    int threads = 0;
#endif
    char *outfname = NULL;
    int c;
    char *infname;
    if(argc < 6) {
        fprintf(stderr, "Too few arguments. See usage.\n");
        print_usage(argv);
        return 1;
    }
#ifdef BMF_THREADS
    while ((c = getopt(argc, argv, "l:o:t:")) > -1) {
#endif
    while ((c = getopt(argc, argv, "l:o:t:")) > -1) {
        switch(c) {
            case 'l': blen = atoi(optarg); break;
#ifdef BMF_THREADS
            case 't': threads = atoi(optarg); break;
#endif
            case 'o': outfname = strdup(optarg); break;
            default: print_opt_err(argv, optarg);
        }
    }
    infname = strdup(argv[optind]);
    if(!blen) {
        fprintf(stderr, "-l option (barcode length) required for bmftools_dmp.\n See usage.\n");
        return 1;
    }
    int retcode = bmftools_dmp_wrapper(infname, outfname, blen);
#if !NDEBUG
    // DEBUG code goes here.
#endif
    free(infname);
    free(outfname);
    if(retcode) {
        fprintf(stderr, "bmftools_dmp_wrapper and bmftools_dmp_core returned a non-zero exit status. (EXIT_FAILURE). Abort!\n");
    }
    return retcode;
}


int bmftools_dmp_wrapper(char *input_path, char *output_path,
                         int blen) {
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
    int retcode = bmftools_dmp_core(seq, out_handle, blen);
    kseq_destroy(seq);
    return retcode;
}


int bmftools_dmp_core(kseq_t *seq, FILE *out_handle, int blen) {
    int l, readlen;
    int *nuc_indices = malloc(2 * sizeof(int));
    char *current_barcode = (char *)calloc(blen + 1, 1);
    l = kseq_read(seq);
    if(l < 0) {
        fprintf(stderr, "Could not open input fastq. Abort!\n");
        exit(1);
    }
    readlen = strlen(seq->seq.s);
    KingFisher_t Holloway = init_kf(readlen);
    KingFisher_t *Hook = &Holloway;
    char *bs_ptr = barcode_mem_view(seq); //Note: NOT null-terminated at the end of the barcode.
    memcpy(current_barcode, bs_ptr, blen);
    pushback_kseq(Hook, seq, nuc_indices); // Initialize Hook with
    while ((l = kseq_read(seq)) >= 0) {
        bs_ptr = barcode_mem_view(seq);
#if !NDEBUG
        if(!bs_ptr) { // bs_ptr is NULL
            destroy_kf(Hook);
            fprintf(stderr, "Malformed fastq comment field - missing the second delimiter. Abort!\n");
            return 1;
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
    return 0;
}
