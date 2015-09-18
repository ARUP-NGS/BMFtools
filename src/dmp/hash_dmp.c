#include "include/kingfisher.h"
#include "dmp_interface.h"  // Contains the KingFisher_t type and KSEQ_INIT
#include "igamc_cephes.c"

extern double igamc(double a, double x);



typedef struct tmpvars {
    char *bs_ptr;
    int blen;
    int readlen;
    int nuc_indices[2];
    int l; // For holding ret value for seq.
} tmpvars_t;


typedef struct HashKing {
    UT_hash_handle hh;
    int64_t id;
    KingFisher_t *value;
} HashKing_t;


extern int64_t get_binnerl(char *barcode, int length); // From binner.h
extern char *barcode_mem_view(kseq_t *seq); // from dmp.h
char ARRG_MAX_TO_NUC(int argmaxret);
double igamc_pvalues(int num_pvalues, double x);
int ARRG_MAX(KingFisher_t *kfp, int index);
int infer_barcode_length(char *bs_ptr);
int pvalue_to_phred(double pvalue);
void destroy_kf(KingFisher_t *kfp);
void dmp_process_write(KingFisher_t *kfp, FILE *handle, int blen, tmpbuffers_t tmp);
void fill_csv_buffer(int readlen, int *arr, char *buffer, char *prefix);
void fill_fa_buffer(KingFisher_t *kfp, int *agrees, char *buffer);
void fill_pv_buffer(KingFisher_t *kfp, int *agrees, char *buffer);
void pushback_kseq(KingFisher_t *kfp, kseq_t *seq, int *nuc_indices, int blen);
KingFisher_t init_kf(int readlen);
void nuc_to_pos(char character, int *nuc_indices);
char test_hp(kseq_t *seq, int threshold);
int get_binner(char *barcode, int length);
void hash_dmp_core(FILE *handle, HashKing_t *hash, tmpvars_t tmp, kseq_t *seq);
int64_t lpow(int64_t base, int64_t exp);
static void pushback_entry(HashKing_t *hash, HashKing_t *tmp_entry, tmpvars_t *tmp, kseq_t *seq);


tmpvars_t init_tmpvars(char *bs_ptr, int blen, int readlen)
{
    tmpvars_t ret = {
            .blen = blen,
            .readlen = readlen,
            .bs_ptr = bs_ptr,
            .nuc_indices = {0, 0},
    };
}



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
    fprintf(stderr, "Now reading from file or handle %s.\n", infname);
    kseq_t *seq = kseq_init(fp);
    fprintf(stderr, "Opened file handles, initiated kseq parser.\n");
    // Initialized kseq
    int l = kseq_read(seq);
    if(l < 0) {
        fprintf(stderr, "Could not open fastq file. Abort mission!\n");
        exit(1);
    }
    fprintf(stderr, "Now getting my barcode pointer.\n");
    char *bs_ptr = barcode_mem_view(seq);
    int blen = infer_barcode_length(bs_ptr);
    int readlen = seq->seq.l;
    fprintf(stderr, "About to start my tmpvars.\n");
    tmpvars_t tmp = init_tmpvars(bs_ptr, blen, readlen);
    HashKing_t *hash = NULL;
    fprintf(stderr, "About to start my hash.\n");
    HashKing_t *current_entry = (HashKing_t *)malloc(sizeof(HashKing_t));
    int64_t bin = get_binnerl(bs_ptr, blen);
    current_entry->id = bin;
    *current_entry->value = init_kf(readlen);
    fprintf(stderr, "Now about to add a value to my hash.\n");
    pushback_kseq(current_entry->value, seq, tmp.nuc_indices, blen);
    fprintf(stderr, "Pushed back kseq.\n");
    HASH_ADD(hh, hash, id, sizeof(int64_t), current_entry);
    fprintf(stderr, "Initiated hash table.\n");
    //pushback_entry(hash, current_entry, &tmp);
    hash_dmp_core(out_handle, hash, tmp, seq);
    fprintf(stderr, "Returned from hash_dmp_core.\n");
    kseq_destroy(seq);
    free(hash);
    fclose(out_handle);
    return 0;
}

static inline void pushback_entry(HashKing_t *hash, HashKing_t *tmp_entry, tmpvars_t *tmp, kseq_t *seq)
{
    tmp->bs_ptr = barcode_mem_view(seq);
    int64_t bin = get_binnerl(tmp->bs_ptr, tmp->blen);
    fprintf(stderr, "Bin for pushing back hash: %i.", bin);
    HASH_FIND(hh, hash, &bin, sizeof(int64_t), tmp_entry);
    if(tmp_entry) {
        pushback_kseq(tmp_entry->value, seq, tmp->nuc_indices, tmp->blen);
    }
    else {
        tmp_entry = (HashKing_t *)malloc(sizeof(tmp_entry));
        tmp_entry->id = bin;
        tmp_entry->value = (KingFisher_t *)malloc(sizeof(KingFisher_t));
        *tmp_entry->value = init_kf(tmp->readlen);
        pushback_kseq(tmp_entry->value, seq, tmp->nuc_indices, tmp->blen);
    }
    return;
}

static inline void pushback_hash(FILE *handle, HashKing_t *hash, tmpvars_t *tmp, kseq_t *seq)
{
    HashKing_t *tmp_entry;
    tmp->bs_ptr = barcode_mem_view(seq);
    int64_t bin = get_binnerl(tmp->bs_ptr, tmp->blen);
    fprintf(stderr, "Bin for pushing back hash: %i.", bin);
    HASH_FIND(hh, hash, &bin, sizeof(int64_t), tmp_entry);
    if(!tmp_entry) { // Entry not found in hash_table.
        fprintf(stderr, "New barcode! %s, %i.", tmp->bs_ptr, bin);
        tmp_entry = (HashKing_t *)malloc(sizeof(HashKing_t));
        tmp_entry->id = bin;
        tmp_entry->value = (KingFisher_t *)malloc(sizeof(KingFisher_t));
        *tmp_entry->value = init_kf(tmp->readlen);
        pushback_kseq(tmp_entry->value, seq, tmp->nuc_indices, tmp->blen);
        HASH_ADD(hh, hash, id, sizeof(int64_t), tmp_entry);
    }
    else {
        pushback_kseq(tmp_entry->value, seq, tmp->nuc_indices, tmp->blen);
    }
    return;
}


void hash_dmp_core(FILE *handle, HashKing_t *hash, tmpvars_t tmp, kseq_t *seq) {
	tmpbuffers_t tmpbuffers;
    int l;
    l = kseq_read(seq);
    HashKing_t *current_entry, *tmp_entry;
    fprintf(stderr, "Now beginning hash_dmp_core.\n");
    int ret;
#if !NDEBUG
    fprintf(stderr, "Seq value. %s", seq->seq.s);
    tmp.bs_ptr = barcode_mem_view(seq);
    fprintf(stderr, "Bs PTR! %s.", tmp.bs_ptr);
    char *first_barcode = (char *)calloc(tmp.blen + 1, sizeof(char));
    first_barcode[tmp.blen] = '\0';
    memcpy(first_barcode, tmp.bs_ptr, tmp.blen);
    fprintf(stderr, "Hey, initialized barcode string is %s.\n", first_barcode);
#endif
    // Delete every between here and "int l" when done.
    while ((l = kseq_read(seq)) >= 0) {
        fprintf(stderr, "Hey, kseq now has seq %s and qual %s.\n", seq->seq.s, seq->qual.s);
        pushback_entry(hash, tmp_entry, &tmp, seq);
    }
    int barcode_counter = 0;

    HASH_ITER(hh, hash, current_entry, tmp_entry) {
        dmp_process_write(current_entry->value, handle, tmp.blen, tmpbuffers);
        free(current_entry->value);
        HASH_DEL(hash, current_entry);
        free(current_entry);
        ++barcode_counter;
    }
    /*
     * Write out results!
     */
    fprintf(stderr, "Number of unique barcodes: %i\n", barcode_counter);
    //dmp_process_write(Hook, handle, Navy.bs_ptr, blen, tmpbuffers)
}
