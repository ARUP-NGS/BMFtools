#include "include/khash.h"
#include "dmp_interface.h"  // Contains the KingFisher_t type and KSEQ_INIT


const int fisher = 137; // I just picked this number?
KHASH_MAP_INIT_INT(fisher, KingFisher_t *) // Sets up khash table for string keys and KingFisher_t * return values.

typedef struct hash_dmp_data {
    khash_t(fisher) *hash_table;
    FILE *out_handle;
    char *out_fname;
    kseq_t *seq;
    khiter_t k;
} hash_dmp_data_t;

// k must be a khiter_t object in scope.
// shorthand way to get the key from hashtable or defVal if not found
#define kh_get_val(kname, hash, key, defVal, k) ({k=kh_get(kname, hash, key);(k!=kh_end(hash)?kh_val(hash,k):defVal);})

// shorthand way to set value in hash with single line command.  Returns value
// returns 0=replaced existing item, 1=bucket empty (new key), 2-adding element previously deleted
#define kh_set(kname, hash, key, val) ({int ret; k = kh_put(kname, hash,key,&ret); kh_value(hash,k) = val; ret;})

extern int get_binner(char *barcode, int length); // From fqmarksplit.h
extern char *barcode_mem_view(kseq_t *seq); // from dmp.h
char ARRG_MAX_TO_NUC(int argmaxret);
double igamc_pvalues(int num_pvalues, double x);
int ARRG_MAX(KingFisher_t *kfp, int index);
int infer_barcode_length(char *bs_ptr);
int pvalue_to_phred(double pvalue);
void destroy_kf(KingFisher_t *kfp);
void dmp_process_write(KingFisher_t *kfp, FILE *handle, int blen);
void fill_csv_buffer(int readlen, int *arr, char *buffer, char *prefix);
void fill_fa_buffer(KingFisher_t *kfp, int *agrees, char *buffer);
void fill_pv_buffer(KingFisher_t *kfp, int *agrees, char *buffer);
void hash_dmp_core(khash_t(fisher) *hash_table, FILE *handle, kseq_t *seq, khiter_t k);
void pushback_hash(khash_t(fisher) *hash_table, kseq_t *seq, char *bs_ptr, int blen, int readlen, int *nuc_indices, khiter_t k);
void pushback_kseq(KingFisher_t *kfp, kseq_t *seq, int *nuc_indices, int blen);

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
    kseq_t *seq = kseq_init(fp);
    fprintf(stderr, "Opened file handles, initiated kseq parser.\n");
    khiter_t k;
    khash_t(fisher) *hash_table = kh_init(fisher);
    fprintf(stderr, "Initiated hash table.\n");
#if !DRY_RUN
    hash_dmp_core(hash_table, out_handle, seq, k);
#endif
    kh_destroy(fisher, hash_table);
    kseq_destroy(seq);
    return 0;
}


void hash_dmp_core(khash_t(fisher) *hash_table, FILE *handle, kseq_t *seq, khiter_t k) {
    int *nuc_indices = malloc(2 * sizeof(int));
    char *bs_ptr = barcode_mem_view(seq);
    int blen = infer_barcode_length(bs_ptr);
    int l = kseq_read(seq);
    int readlen = strlen(seq->seq.s);
    fprintf(stderr, "read length (inferred): %i\n", readlen);
    if(l < 0) {
        fprintf(stderr, "Could not open fastq file. Abort mission!\n");
        exit(1);
    }
    pushback_hash(hash_table, seq, bs_ptr, blen, readlen, nuc_indices, k);
    while ((l = kseq_read(seq)) >= 0) {
        bs_ptr = barcode_mem_view(seq);
        pushback_hash(hash_table, seq, bs_ptr, blen, readlen, nuc_indices, k);
    }
    KingFisher_t *Hook; // KingFisher pointer through which to iterate.
    kh_foreach_value(hash_table, Hook, {
            dmp_process_write(Hook, handle, blen);
            destroy_kf(Hook);
        });
    //dmp_process_write(Hook, handle, bs_ptr, blen)
    free(nuc_indices);
}

//int parallel_hash_dmp(char )

void pushback_hash(khash_t(fisher) *hash_table, kseq_t *seq, char *bs_ptr, int blen, int readlen, int *nuc_indices, khiter_t k) {
    KingFisher_t *Hook = kh_get_val(fisher, hash_table, get_binner(bs_ptr, blen), NULL, k);
    if(!Hook) { // Hash table was empty - start a new Holloway!
        INIT_KF(Holloway, readlen)
        Hook = &Holloway;
    }
    pushback_kseq(Hook, seq, nuc_indices, blen);
}
