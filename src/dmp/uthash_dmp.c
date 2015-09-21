#include "include/kingfisher.h"
#include "dmp_interface.h"  // Contains the KingFisher_t type and KSEQ_INIT
#include "uthash_dmp.h"
#include "igamc_cephes.c"

extern double igamc(double a, double x);



typedef struct tmpvars {
    char *bs_ptr;
    int blen;
    int readlen;
    int nuc_indices[2];
    char key[MAX_BARCODE_LENGTH + 1];
    int l; // For holding ret value for seq.
} tmpvars_t;


typedef struct HashKing {
    UT_hash_handle hh;
    char id[MAX_BARCODE_LENGTH + 1];
    KingFisher_t *value;
} HashKing_t;


extern uint64_t get_binnerul(char *barcode, int length); // From binner.h
extern int64_t get_binnerl(char *barcode, int length); // From binner.h
extern char *barcode_mem_view(kseq_t *seq); // from dmp.h
int64_t lpow(int64_t base, int64_t exp);
char ARRG_MAX_TO_NUC(int argmaxret);
double igamc_pvalues(int num_pvalues, double x);
int ARRG_MAX(KingFisher_t *kfp, int index);
int infer_barcode_length(char *bs_ptr);
int pvalue_to_phred(double pvalue);
void destroy_kf(KingFisher_t *kfp);
void dmp_process_write(KingFisher_t *kfp, FILE *handle, int blen, tmpbuffers_t *tmp);
void fill_csv_buffer(int readlen, int *arr, char *buffer, char *prefix, char typecode);
void fill_fa_buffer(KingFisher_t *kfp, int *agrees, char *buffer);
void fill_pv_buffer(KingFisher_t *kfp, int *agrees, char *buffer);
void pushback_kseq(KingFisher_t *kfp, kseq_t *seq, int *nuc_indices, int blen);
KingFisher_t init_kf(int readlen);
void nuc_to_pos(char character, int *nuc_indices);
char test_hp(kseq_t *seq, int threshold);
int get_binner(char *barcode, int length);
void hash_dmp_core(FILE *handle, HashKing_t *hash, tmpvars_t tmp, kseq_t *seq);
uint64_t ulpow(uint64_t base, uint64_t exp);
void cp_view2buf(char *view, char *buf);
void omgz_core(FILE *handle, HashKing_t *hash, tmpvars_t *tmp, kseq_t *seq);
void cp_view2buf(char *view, char *buf);


tmpvars_t init_tmpvars(char *bs_ptr, int blen, int readlen)
{
    tmpvars_t ret = {
            .blen = blen,
            .readlen = readlen,
            .bs_ptr = bs_ptr,
            .nuc_indices = {0, 0}
    };
    return ret;
}



void print_usage(char *argv[]) {
    fprintf(stderr, "Usage: %s -o <output_filename> <input_filename>.\n", argv[0]);
}

void print_opt_err(char *argv[], char *optarg) {
    fprintf(stderr, "Invalid argument %s. See usage.\n", optarg);
    print_usage(argv);
    exit(1);
}


inline void cp_view2buf(char *view, char *buf)
{
    int blen = 0;
    while(view[blen] != '\0' && view[blen] != '|') {
        ++blen;
    }
    memcpy(buf, view, blen);
    buf[blen] = '\0';
    return;
}

int main(int argc, char *argv[]){
    char *outfname = NULL;
    char *infname = NULL;
    int c;
    while ((c = getopt(argc, argv, "ho:")) > -1) {
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
    fprintf(stderr, "Barcode length: %i.\n", blen);
#if !NDEBUG
    char first_bc[MAX_BARCODE_LENGTH + 1];
    cp_view2buf(bs_ptr, first_bc);
#endif
    tmpvars_t tmp = init_tmpvars(bs_ptr, blen, seq->seq.l);
    fprintf(stderr, "tmpvars: bs_ptr %s, blen %i, readlen %i.\n", bs_ptr, tmp.blen, tmp.readlen);
    HashKing_t *hash = NULL;
    HashKing_t *current_entry = (HashKing_t *)calloc(1, sizeof(HashKing_t));
    HashKing_t *tmp_hk = current_entry; // Save the pointer location for later comparison.
    cp_view2buf(bs_ptr, current_entry->id);
    fprintf(stderr, "About to start my hash with key %s and readlen %i.\n", current_entry->id, tmp.readlen);
    current_entry->value = init_kfp(tmp.readlen);
    HASH_ADD_STR(hash, id, current_entry);
    fprintf(stderr, "Now about to add a value to my kf, whose pointer is at %p.\n", current_entry->value);
    pushback_kseq(current_entry->value, seq, tmp.nuc_indices, tmp.blen);
    fprintf(stderr, "Pushed back kseq.\n");
    fprintf(stderr, "Initiated hash table.\n");
    /*
     * Start of testing block
     */

    tmpbuffers_t tmpbuffers;
    omgz_core(out_handle, hash, &tmp, seq);
#if !NDEBUG
    HASH_FIND_STR(hash, first_bc, current_entry);
    fprintf(stderr, "Fetched entry %p which should have the same value as %p.\n", current_entry, tmp_hk);
#endif
    fprintf(stderr, "Loaded all fastq records into memory for meta-analysis. Now writing out to file!\n");
    HASH_ITER(hh, hash, current_entry, tmp_hk) {
#if !NDEBUG
        fprintf(stderr, "Now writing collapsed family with barcode %s.\n", current_entry->id);
#endif
        dmp_process_write(current_entry->value, out_handle, tmp.blen, &tmpbuffers);
        destroy_kf(current_entry->value);
        free(current_entry->value);
        current_entry->value = NULL;
        HASH_DEL(hash, current_entry);
        free(current_entry);
    }
    //hash_dmp_core(out_handle, hash, tmp, seq);
    //fprintf(stderr, "Returned from hash_dmp_core.\n");
    kseq_destroy(seq);
    free(hash);
    fclose(out_handle);
    free(outfname);
    free(infname);
    return 0;
}

void omgz_core(FILE *handle, HashKing_t *hash, tmpvars_t *tmp, kseq_t *seq)
{
    HashKing_t *current_entry;
    int l;
    while((l = kseq_read(seq)) >= 0) {
        tmp->bs_ptr = barcode_mem_view(seq);
        cp_view2buf(tmp->bs_ptr, tmp->key);
        HASH_FIND_STR(hash, tmp->key, current_entry);
        if(!current_entry) {
#if !NDEBUG
            fprintf(stderr, "Good, this entry shouldn't have been found.\n");
#endif
            current_entry = (HashKing_t *)malloc(sizeof(HashKing_t));
            current_entry->value = init_kfp(tmp->readlen);
            cp_view2buf(tmp->bs_ptr, current_entry->id);
            pushback_kseq(current_entry->value, seq, tmp->nuc_indices, tmp->blen);
            HASH_ADD_STR(hash, id, current_entry);
        }
        else {
#if !NDEBUG
            fprintf(stderr, "Barcode found already (id: %s.).\n", current_entry->id);
#endif
            pushback_kseq(current_entry->value, seq, tmp->nuc_indices, tmp->blen);
        }
    }

}

