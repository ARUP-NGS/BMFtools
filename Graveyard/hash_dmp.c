#include "include/kingfisher.h"
#include "dmp_interface.h"  // Contains the KingFisher_t type and KSEQ_INIT

extern int64_t get_binnerl(char *barcode, int length); // From binner.h
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
void hash_dmp_core(outpost_t Navy, FILE *handle);
void pushback_kseq(KingFisher_t *kfp, kseq_t *seq, int *nuc_indices, int blen);
KingFisher_t init_kf(int readlen);
void nuc_to_pos(char character, int *nuc_indices);
char test_hp(kseq_t *seq, int threshold);
static void pushback_hash(outpost_t Navy);
int get_binner(char *barcode, int length);
int64_t lpow(int64_t base, int64_t exp);



void print_usage(char *argv[]) {
    fprintf(stderr, "Usage: %s -o <output_filename> <input_filename>.\n", argv[0]);
}

void print_opt_err(char *argv[], char *optarg) {
    fprintf(stderr, "Invalid argument %c. See usage.\n", optarg);
    print_usage(argv);
    exit(1);
}


static void pushback_hash(outpost_t Navy)
{
    Navy.bs_ptr = barcode_mem_view(Navy.seq);
    int64_t bin = get_binnerl(Navy.bs_ptr, Navy.blen);
    fprintf(stderr, "Bin for pushing back hash: %i.", bin);
    Navy.k=kh_get(fisher, Navy.hash,
                  get_binner(Navy.bs_ptr, Navy.blen));
    if(Navy.k==kh_end(Navy.hash)) {
    	fprintf(stderr, "New barcode! %s, %i.", Navy.bs_ptr, bin);
        set_kf(Navy.seq->seq.l, *kh_val(Navy.hash, Navy.k));
        pushback_kseq(kh_val(Navy.hash, Navy.k), Navy.seq, Navy.nuc_indices, Navy.blen);
    }
    else {
        pushback_kseq(kh_val(Navy.hash, Navy.k), Navy.seq, Navy.nuc_indices, Navy.blen);
    }
    return;
}

int main(int argc, char *argv[]){
    char *outfname = NULL;
    char *infname = NULL;
    int c;
    while ((c = getopt(argc, argv, "h:o:")) > -1) {
        switch(c) {
            case 'o': outfname = strdup(optarg);break;
            case 'h': print_usage(argv); return 0;
            default: print_opt_err(argv, optarg);
        }
    }
    if(argc < 4) {
        fprintf(stderr, "Required arguments missing. See usage.\n");
        print_usage(argv);
        exit(1);
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
    int ret = -1;
    fprintf(stderr, "Trying to start hash table.\n");
    k = kh_put(fisher, hash, 133, &ret);
    outpost_t Navy = {
            .hash = hash,
            .seq = seq,
            .bs_ptr = bs_ptr,
            .blen = blen,
            .k = k,
            .nuc_indices = nuc_indices
    };
    hash_dmp_core(Navy, out_handle);
    fprintf(stderr, "Returned from hash_dmp_core.\n");
    kh_destroy(fisher, hash);
    kseq_destroy(seq);
    fclose(out_handle);
    return 0;
}


void hash_dmp_core(outpost_t Navy, FILE *handle) {
	khint_t k;
    int ret = 0;
    int64_t bin = get_binnerl(Navy.bs_ptr, Navy.blen);
    fprintf(stderr, "Now beginning hash_dmp_core.\n");
    k = kh_put(fisher, Navy.hash, bin, &ret);
    fprintf(stderr, "New KingFisher_t! Barcode: %s. Bin: %i.", kh_val(Navy.hash, k)->barcode, bin);
    if(!Navy.bs_ptr) {
        fprintf(stderr, "Locating barcode failed. seq's comment: %s.\n", Navy.seq->comment.s);
        exit(1);
    }
    k = kh_put(fisher, Navy.hash, bin, &ret);
    fprintf(stderr, "New bin %s.\n", k);
    *kh_val(Navy.hash, k) = init_kf(Navy.seq->seq.l);
    fprintf(stderr, "New KingFisher_t! Barcode: %s. Bin: %i.", kh_val(Navy.hash, k)->barcode, bin);
    pushback_kseq(kh_val(Navy.hash, k), Navy.seq, Navy.nuc_indices, Navy.blen);
    /*
    char *first_barcode = (char *)malloc((Navy.blen + 1) * sizeof(char));
    memcpy(first_barcode, Navy.bs_ptr, Navy.blen);
    first_barcode[Navy.blen] = '\0';
    */
    //pushback_hash(hash, Navy.seq, Navy.bs_ptr, Navy.blen, readlen, nuc_indices, k);
    // Delete every between here and "int l" when done.
    int l;
    while ((l = kseq_read(Navy.seq)) >= 0) {
        pushback_hash(Navy);
    }
    int barcode_counter = 0;
    khint_t _i;
    for(_i = kh_begin(Navy.hash); _i != kh_end(Navy.hash); ++_i) {
        if(!kh_exist(Navy.hash, _i)) continue;
        fprintf(stderr, "Barcode for this KingFisher_t: %s.\n", kh_val(Navy.hash, _i)->barcode);
        fprintf(stderr, "Just got the Hook for the KingFisher_t object.\n");
        fprintf(stderr, "Now let's just destroy Hook.\n");
    }
    for(_i = kh_begin(Navy.hash); _i != kh_end(Navy.hash); ++_i) {
        if(!kh_exist(Navy.hash, _i)) continue;
        dmp_process_write(kh_val(Navy.hash, _i), handle, Navy.blen);
        destroy_kf(kh_val(Navy.hash, _i));
        barcode_counter++;
    }
    fprintf(stderr, "Number of unique barcodes: %i\n", barcode_counter);
    //dmp_process_write(Hook, handle, Navy.bs_ptr, blen)
    free(Navy.nuc_indices);
}
