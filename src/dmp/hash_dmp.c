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
void pushback_kseq(KingFisher_t *kfp, kseq_t *seq, int *nuc_indices, int blen);

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
    int blen; // Barcode length;
    int l; // The kseq iterator integer
    int ret; // The khash return integer
    char *bs_ptr; // The pointer to the Barcode Sequence
    int *nuc_indices; // Temporary variable for storing indices for accessing the KingFisher arrays
    uint64_t bin; // bin to place the index in
} Outpost_t;

Outpost_t *INIT_OUTPOST(FILE *in_handle, FILE *out_handle) {
    gzFile fp = gzdopen(fileno(in_handle), "r");
    kseq_t *seq = kseq_init(fp);
    int l = kseq_read(seq);
    if(l < 0) {
        fprintf(stderr, "Couldn't open input file handle for reading. Abort!\n");
        exit(1);
    }
    int readlen = strlen(seq->seq.s);
    char *bs_ptr = barcode_mem_view(seq);
    int blen = infer_barcode_length(bs_ptr);
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
        .out_handle = out_handle,
        .l = 0,
        .ret = 0,
        .bs_ptr = NULL,
        .blen = blen,
        .nuc_indices = (int *)malloc(2 * sizeof(int)),
        .bin = 0
    };
    Outpost_t *ret_ptr = &ret;
    return ret_ptr;
}


void RESIZE_OUTPOST(Outpost_t *outpost) {
    uint64_t tmpmax = outpost->max;
    outpost->max += EXPEDITION_INC_SIZE;

    outpost->Expedition = (KingFisher_t **)realloc(outpost->Expedition, (outpost->max + EXPEDITION_INC_SIZE) * sizeof(KingFisher_t *));
    if(!outpost->Expedition) { // Realloc failed!
        fprintf(stderr, "Reallocating memory for Outpost to final size %i failed. Abort!\n", outpost->max);
        exit(1);
    }
    for(int i = outpost->max - EXPEDITION_INC_SIZE; i < outpost->max; ++i) {
        outpost->Expedition[i] = init_kfp(outpost->Expedition[0]->readlen);
    }
    return;
}


/*
 * @returns: 1 upon success, 0 upon end of file.
 */
void PUSHBACK_OUTPOST(Outpost_t *outpost) {
    outpost->bs_ptr = barcode_mem_view(outpost->seq);
    outpost->bin = get_binnerl(outpost->bs_ptr, outpost->blen);
    outpost->k = kh_get(fisher, outpost->hash, outpost->bin);
    if(outpost->k == kh_end(outpost->hash)) {
        outpost->k = kh_put(fisher, outpost->hash, outpost->bin, &(outpost->ret));
        kh_value(outpost->hash, outpost->k) = outpost->length;
        pushback_kseq(outpost->Expedition[outpost->length - 1], outpost->seq,
                      outpost->nuc_indices, outpost->blen);
        outpost->length++;
        if(outpost->length == outpost->max) {
            RESIZE_OUTPOST(outpost);
        }
    }
    else {
        pushback_kseq(outpost->Expedition[kh_value(outpost->hash, outpost->k)],
                      outpost->seq,
                      outpost->nuc_indices, outpost->blen);
    }
    return;
}

/*
 * Cleans up temporary structures. Does not close either file handle.
 * Expedition's entries also need to have been cleaned up before this is called!
 */

void destroy_outpost(Outpost_t *outpost) {
    free(outpost->bs_ptr);
    free(outpost->nuc_indices);
    kseq_destroy(outpost->seq);
    kh_destroy(fisher, outpost->hash);
    free(outpost->Expedition);
    return;
}


void write_expedition_results(Outpost_t *outpost) {
    for(int i = 0; i < outpost->length; i++) {
        dmp_process_write(outpost->Expedition[i], outpost->out_handle, outpost->blen);
        destroy_kf(outpost->Expedition[i]);
    }
    for(int i = outpost->length; i < outpost->max; ++i) {
        destroy_kf(outpost->Expedition[i]); // Get rid of the pre-allocated but unused KingFisher_t objects.
    }
}

void outpost_dmp_core(FILE *in_handle, FILE *out_handle) {
    Outpost_t *outpost = INIT_OUTPOST(in_handle, out_handle);
    PUSHBACK_OUTPOST(outpost); // Run this just once before the loop, as I've iterated through it once to get read length and barcod length.
    while((outpost->l = kseq_read(outpost->seq)) >= 0) {
        PUSHBACK_OUTPOST(outpost);
    }
    write_expedition_results(outpost);
    destroy_outpost(outpost);
}



//TODO: - to fix this attempted khash use.
// 7. ???
// 8. PROFIT

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
    outpost_dmp_core(in_handle, out_handle);
    fprintf(stderr, "Returned from outpost_dmp_core.\n");
    fclose(out_handle);
    fclose(in_handle);
    return 0;
}
