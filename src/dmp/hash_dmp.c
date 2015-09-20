#include "include/kingfisher.h"
#include "dmp_interface.h"  // Contains the KingFisher_t type and KSEQ_INIT
#include "igamc_cephes.c"

extern int64_t get_binnerl(char *barcode, int length); // From binner.h
extern char *barcode_mem_view(kseq_t *seq); // from dmp.h
char ARRG_MAX_TO_NUC(int argmaxret);
double igamc_pvalues(int num_pvalues, double x);
int ARRG_MAX(KingFisher_t *kfp, int index);
int infer_barcode_length(char *bs_ptr);
int pvalue_to_phred(double pvalue);
void destroy_kf(KingFisher_t *kfp);
void dmp_process_write(KingFisher_t *kfp, FILE *handle, int blen, tmpbuffers_t *tmp);
void fill_csv_buffer(int readlen, int *arr, char *buffer, char *prefix);
void fill_fa_buffer(KingFisher_t *kfp, int *agrees, char *buffer);
void fill_pv_buffer(KingFisher_t *kfp, int *agrees, char *buffer);
void hash_dmp_core(outpost_t *Navy, FILE *handle);
void pushback_kseq(KingFisher_t *kfp, kseq_t *seq, int *nuc_indices, int blen);
KingFisher_t init_kf(int readlen);
void nuc_to_pos(char character, int *nuc_indices);
char test_hp(kseq_t *seq, int threshold);
void pushback_hash(outpost_t *Navy);
int get_binner(char *barcode, int length);
uint64_t get_binnerul(char *barcode, int length);
uint64_t ulpow(uint64_t base, uint64_t exp);
int64_t lpow(int64_t base, int64_t exp);
//extern double igamc(double x, double y);



void print_usage(char *argv[]) {
    fprintf(stderr, "Usage: %s -o <output_filename> <input_filename>.\n", argv[0]);
}

void print_opt_err(char *argv[], char *optarg) {
    fprintf(stderr, "Invalid argument %s. See usage.\n", optarg);
    print_usage(argv);
    exit(1);
}

int main(int argc, char *argv[]){
    char *outfname = NULL;
    char *infname = NULL;
    int c;
    if(argc < 4) {
        fprintf(stderr, "Required arguments missing. See usage.\n");
        print_usage(argv);
        exit(1);
    }
    while ((c = getopt(argc, argv, "h:o:")) > -1) {
        switch(c) {
            case 'o': outfname = strdup(optarg);break;
            case 'h': print_usage(argv); return 0;
            default: print_opt_err(argv, optarg);
        }
    }
    FILE *in_handle;
    FILE *out_handle;
    infname = strdup(argv[optind]);
    if(infname[0] == '-' || !infname) in_handle = stdin;
    else {
        in_handle = fopen(infname, "r");
    }
    if(!outfname) out_handle = stdout;
    else {
        out_handle = fopen(outfname, "w");
    }
    gzFile fp = gzdopen(fileno(in_handle), "r");
    fprintf(stderr, "Now reading from file or handle %s", infname);
    kseq_t *seq = kseq_init(fp);
    fprintf(stderr, "Opened file handles, initiated kseq parser.\n");
    // Initialized kseq
    int l = kseq_read(seq);
    if(l < 0) {
        fprintf(stderr, "Could not open fastq file. Abort mission!\n");
        exit(1);
    }
    char *bs_ptr = barcode_mem_view(seq);
    int blen = infer_barcode_length(bs_ptr);
    int readlen = seq->seq.l;
    int *nuc_indices = (int *)calloc(2, sizeof(int));
    khiter_t k;
    khash_t(fisher) *hash = kh_init(fisher);
    fprintf(stderr, "Initiated hash table.\n");
    outpost_t Navy = {
            .hash = hash,
            .seq = seq,
            .bs_ptr = bs_ptr,
            .blen = blen,
            .k = k,
            .nuc_indices = nuc_indices,
            .key = 0,
            .ret = 0
    };
    hash_dmp_core(&Navy, out_handle);
    fprintf(stderr, "Returned from hash_dmp_core.\n");
    kh_destroy(fisher, hash);
    kseq_destroy(seq);
    fclose(out_handle);
    return 0;
}


void hash_dmp_core(outpost_t *Navy, FILE *handle) {
    fprintf(stderr, "Now beginning hash_dmp_core.\n");
    if(!Navy->bs_ptr) {
        fprintf(stderr, "Locating barcode failed. seq's comment: %s.\n", Navy->seq->comment.s);
        exit(1);
    }
    KingFisher_t *Holloway = init_kfp(Navy->seq->seq.l);
    //fprintf(stderr, "Holloway's current length: %i. Barcode: %s. Pointer to Holloway: %p.\n", Holloway.length, Holloway.barcode, &Holloway);
    pushback_kseq(Holloway, Navy->seq, Navy->nuc_indices, Navy->blen);
    int ret;
    Navy->key = get_binnerul(Navy->bs_ptr, Navy->blen);
    Navy->k = kh_put(fisher, Navy->hash, Navy->key, &ret);
    kh_value(Navy->hash, Navy->k) = Holloway;
    /*
    char *first_barcode = (char *)malloc((Navy->blen + 1) * sizeof(char));
    memcpy(first_barcode, Navy->bs_ptr, Navy->blen);
    first_barcode[Navy->blen] = '\0';
    */
    //pushback_hash(hash, Navy->seq, Navy->bs_ptr, Navy->blen, readlen, nuc_indices, k);
    // Delete every between here and "int l" when done.
#if !NDEBUG
    fprintf(stderr, "Holloway's current length: %i. Barcode: %s. Pointer to struct: %p.\n", Holloway->length, Holloway->barcode, Holloway);
#endif
    int l;
    while ((l = kseq_read(Navy->seq)) >= 0) {
        pushback_hash(Navy);
    }
    int barcode_counter = 0;
    khint_t _i;
    tmpbuffers_t tmp;
    for(_i = kh_begin(Navy->hash); _i != kh_end(Navy->hash); ++_i) {
        if(!kh_exist(Navy->hash, _i)) continue;
        fprintf(stderr, "Just got the Hook for the KingFisher_t object.\n");
        dmp_process_write(kh_val(Navy->hash, _i), handle, Navy->blen, &tmp);
        fprintf(stderr, "Now let's just destroy Hook.\n");
        //destroy_kf(Hook);
        barcode_counter++;
    }
    fprintf(stderr, "Number of unique barcodes: %i\n", barcode_counter);
    //dmp_process_write(Hook, handle, Navy->bs_ptr, blen)
    free(Navy->nuc_indices);
}
