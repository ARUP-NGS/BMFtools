#include "include/khash.h"
#include "dmp_interface.h"  // Contains the KingFisher_t type and KSEQ_INIT

#define EXPEDITION_INC_SIZE 65536


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
KingFisher_t *init_kfp(int readlen);
const int fisher = 137; // I just picked a number
KHASH_MAP_INIT_INT64(fisher, uint64_t) // Set up khash with 64-bit integer keys and uint64_t payload

typedef struct Outpost {
    KingFisher_t **Expedition; // The KingFishers.
    uint64_t max; // The maximum number of KingFisher_t * objects that can be held before resizing.
    uint64_t length; // The number of KingFisher_t * objects currently filled
    khash_t(fisher) *hash; // My hash table with int64 keys and uint64_t payloads.
    khiter_t k; // The return variable for kh_get
    kseq_t *seq; // The fastq handle
    FILE *out_handle; // The output handle
} Outpost_t;

Outpost_t *INIT_OUTPOST(FILE *in_handle, FILE *out_handle, int readlen) {
    gzFile fp = gzdopen(fileno(in_handle), "r");
    kseq_t *seq = kseq_init(fp);
    khash_t(fisher) *hash = kh_init(fisher);
    KingFisher_t **Expedition = (KingFisher_t **)malloc(EXPEDITION_INC_SIZE * sizeof(KingFisher_t *));
    for(int i = 0; i < EXPEDITION_INC_SIZE; ++i) {
        Expedition[0] = init_kfp(readlen);
    }
    Outpost_t ret = {
        .Expedition = Expedition,
        .max = EXPEDITION_INC_SIZE,
        .hash = hash,
        .length = 0,
        .k = 0,
        .seq = seq,
        .out_handle = out_handle
    };
    Outpost_t *ret_ptr = &ret;
    return ret_ptr;
}

void RESIZE_OUTPOST(Outpost_t *outpost) {
    outpost->max += EXPEDITION_INC_SIZE;
    outpost->Expedition = (KingFisher_t **)realloc(outpost->Expedition, outpost->max * sizeof(KingFisher_t *));
    if(!outpost->Expedition) { // Realloc failed!
        fprintf(stderr, "Reallocating memory for Outpost to final size %i failed. Abort!\n", outpost->max);
        exit(1);
    }
    for(int i = outpost->max - EXPEDITION_INC_SIZE; i < outpost->max; ++i) {
        outpost->Expedition[i] = init_kfp(readlen);
    }
    return;
}

/*
int outpost_dmp_core(Outpost_t *outpost, )
*/

//TODO: - to fix this attempted khash use.
// 1. Write a resize.
// 2. Write a pushback.
// 3. Write a macro each for doing the full put/get in one command.
// 4. Use the old pushback code.
// 5. ???
// 6. PROFIT

void hash_dmp_core(khash_t(fisher) *hash, FILE *handle, kseq_t *seq, khiter_t k);
void pushback_hash(khash_t(fisher) *hash, kseq_t *seq, char *bs_ptr, int blen, int readlen, int *nuc_indices, khiter_t k);
void pushback_kseq(KingFisher_t *kfp, kseq_t *seq, int *nuc_indices, int blen);

void print_usage(char *argv[]) {
    fprintf(stderr, "Usage: %s -o <output_filename> <input_filename>.\n"
                    "Set input_filename to \"-\" to select to stdout.\n"
                    "Leave output_filename unset to select stdin.\n", argv[0]);
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
    if(argc < 2) {
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
    fprintf(stderr, "Now reading from file or handle %s\n", infname);
    kseq_t *seq = kseq_init(fp);
    khiter_t k;
    khash_t(fisher) *hash = kh_init(fisher);
    KingFisher_t **Expedition = (KingFisher_t **)malloc(EXPEDITION_INC_SIZE * sizeof(KingFisher_t *));
    // Initialized kseq
    int l = kseq_read(seq);
    if(l < 0) {
        fprintf(stderr, "Could not open fastq file. Abort mission!\n");
        return 1;
    }
    hash_dmp_core(hash, out_handle, seq, k);
    fprintf(stderr, "Returned from hash_dmp_core.\n");
    kh_destroy(fisher, hash);
    kseq_destroy(seq);
    fclose(out_handle);
    return 0;
}


void hash_dmp_core(khash_t(fisher) *hash, FILE *handle, kseq_t *seq, khiter_t k) {
    fprintf(stderr, "Now beginning hash_dmp_core.\n");
    int *nuc_indices = malloc(2 * sizeof(int));
    char *bs_ptr = barcode_mem_view(seq);
    if(!bs_ptr) {
        fprintf(stderr, "Locating barcode failed. seq's comment: %s.\n", seq->comment.s);
        exit(1);
    }
    int blen = infer_barcode_length(bs_ptr);
    fprintf(stderr, "Now getting barcode length.\n");
    int readlen = strlen(seq->seq.s);
    fprintf(stderr, "read length (inferred): %i\n", readlen);
    KingFisher_t Holloway = init_kf(readlen);
    //fprintf(stderr, "Holloway's current length: %i. Barcode: %s. Pointer to Holloway: %p.\n", Holloway.length, Holloway.barcode, &Holloway);
    pushback_kseq(&Holloway, seq, nuc_indices, blen);
    int ret;
    uint64_t bin = get_binnerl(bs_ptr, blen);
    k = kh_put(fisher, hash, get_binnerl(bs_ptr, blen), &ret);
    /*
    kh_value(hash, k) = &Holloway;
    char *first_barcode = (char *)malloc((blen + 1) * sizeof(char));
    memcpy(first_barcode, bs_ptr, blen);
    first_barcode[blen] = '\0';
    //pushback_hash(hash, seq, bs_ptr, blen, readlen, nuc_indices, k);
    // Delete every between here and "int l" when done.
    fprintf(stderr, "Holloway's current length: %i. Barcode: %s. Pointer to Holloway: %p.\n", Holloway.length, Holloway.barcode, &Holloway);
    int l;
    while ((l = kseq_read(seq)) >= 0) {
        bs_ptr = barcode_mem_view(seq);
        pushback_hash(hash, seq, bs_ptr, blen, readlen, nuc_indices, k);
    }
    k = kh_get(fisher, hash, bin);
    KingFisher_t *omgz = kh_value(hash, k);
    fprintf(stderr, "First barcode's info %i. Barcode: %s. Pointer to Holloway: %p.\n", omgz->length, omgz->barcode, omgz);
    */
    int barcode_counter = 0;
    khint_t _i;
    /*
    for(_i = kh_begin(hash); _i != kh_end(hash); ++_i) {
        if(!kh_exist(hash, _i)) continue;
        fprintf(stderr, "Just got the Hook for the KingFisher_t object.\n");
        dmp_process_write(kh_val(hash, _i), handle, blen);
        fprintf(stderr, "Now let's just destroy Hook.\n");
        //destroy_kf(Hook);
        barcode_counter++;
    }
    */
    fprintf(stderr, "Number of unique barcode: %i\n", barcode_counter);
    //dmp_process_write(Hook, handle, bs_ptr, blen)
    free(nuc_indices);
}

//int parallel_hash_dmp(char )

void pushback_hash(khash_t(fisher) *hash, kseq_t *seq, char *bs_ptr, int blen, int readlen, int *nuc_indices, khiter_t k) {
    //KingFisher_t *Hook = kh_get_val(fisher, hash, get_binner(bs_ptr, blen), NULL, k);
    /*
    fprintf(stderr, "Read length: %i.\n", readlen);
    int ret;
    k = kh_get(fisher, hash, get_binner(bs_ptr, blen));
    if(k == kh_end(hash)) {
        KingFisher_t Holloway = init_kf(readlen);
        k = kh_put(fisher, hash, get_binner(bs_ptr, blen), &ret);
        pushback_kseq(&Holloway, seq, nuc_indices, blen);
        kh_value(hash, k) = &Holloway;
        fprintf(stderr, "Read length for the hash table after being put in: %i\n", Holloway.readlen);
        fprintf(stderr, "Length, barcode after being put in: %i, %s\n", Holloway.length, Holloway.barcode);
        fprintf(stderr, "Value of pointer: %p\n", &Holloway);
    }
    else {
        fprintf(stderr, "Read length of accessed value: %i\n", kh_value(hash, k)->readlen);
        fprintf(stderr, "Length of accessed value: %i\n", kh_value(hash, k)->length);
        fprintf(stderr, "Value of pointer: %p\n", kh_value(hash, k));
        pushback_kseq(kh_value(hash, k), seq, nuc_indices, blen);
    }
    */
    return;
}
