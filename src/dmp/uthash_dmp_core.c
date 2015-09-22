#include "include/kingfisher.h"
#include "dmp_interface.h"  // Contains the KingFisher_t type and KSEQ_INIT
#include "uthash_dmp_core.h"
#include "igamc_cephes.c"

extern double igamc(double a, double x);
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
void omgz_core(char *infname, char *outfname);
void cp_view2buf(char *view, char *buf);


void tmpvars_destroy(tmpvars_t *tmp)
{
    free(tmp->buffers);
    free(tmp);
    tmp = NULL;
    return;
}


tmpvars_t *init_tmpvars_p(char *bs_ptr, int blen, int readlen)
{
    tmpvars_t *ret = (tmpvars_t *)malloc(sizeof(tmpvars_t));
    ret->blen = blen;
    ret->readlen = readlen;
    ret->bs_ptr = bs_ptr;
    ret->nuc_indices[0] = ret->nuc_indices[1] = 0;
    ret->buffers = (tmpbuffers_t *)malloc(sizeof(tmpbuffers_t));
    return ret;
}


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



void print_hashdmp_usage(char *argv[]) {
    fprintf(stderr, "Usage: %s -o <output_filename> <input_filename>.\n", argv[0]);
}

void print_hashdmp_opt_err(char *argv[], char *optarg) {
    fprintf(stderr, "Invalid argument %s. See usage.\n", optarg);
    print_hashdmp_usage(argv);
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


int uthash_dmp_main(int argc, char *argv[])
{
    char *outfname = NULL;
    char *infname = NULL;
    int c;
    while ((c = getopt(argc, argv, "ho:")) > -1) {
        switch(c) {
            case 'o': outfname = strdup(optarg);break;
            case 'h': print_hashdmp_usage(argv); return 0;
            default: print_hashdmp_opt_err(argv, optarg);
        }
    }
    if(argc < 4) {
        fprintf(stderr, "Required arguments missing. See usage.\n");
        print_hashdmp_usage(argv);
        exit(1);
    }
    infname = strdup(argv[optind]);

    omgz_core(infname, outfname);
    free(outfname);
    free(infname);
    return 0;
}

void omgz_core(char *infname, char *outfname)
{
    FILE *in_handle;
    FILE *out_handle;
    if(infname[0] == '-' || !infname) in_handle = stdin;
    else {
        in_handle = fopen(infname, "r");
    }
    fprintf(stderr, "Now reading from file or handle %s.\n", strcmp(infname, "-") == 0 ? "stdin": infname);
    if(!outfname) out_handle = stdout;
    else {
        out_handle = fopen(outfname, "w");
    }
    gzFile fp = gzdopen(fileno(in_handle), "r");
    kseq_t *seq = kseq_init(fp);
    fprintf(stderr, "Opened file handles, initiated kseq parser.\n");
    // Initialized kseq
    int l = kseq_read(seq);
    if(l < 0) {
        fprintf(stderr, "Could not open fastq file (%s). Abort mission!\n",
                strcmp(infname, "-") == 0 ? "stdin": infname);
        exit(1);
    }
    char *bs_ptr = barcode_mem_view(seq);
    int blen = infer_barcode_length(bs_ptr);
    fprintf(stderr, "Barcode length: %i.\n", blen);
    tmpvars_t *tmp = init_tmpvars_p(bs_ptr, blen, seq->seq.l);
    // Start hash table
    HashKing_t *hash = NULL;
    HashKing_t *current_entry = (HashKing_t *)malloc(sizeof(HashKing_t));
    HashKing_t *tmp_hk = current_entry; // Save the pointer location for later comparison.
    cp_view2buf(bs_ptr, current_entry->id);
    fprintf(stderr, "About to start my hash with key %s and readlen %i.\n", current_entry->id, tmp->readlen);
    current_entry->value = init_kfp(tmp->readlen);
    HASH_ADD_STR(hash, id, current_entry);
    pushback_kseq(current_entry->value, seq, tmp->nuc_indices, tmp->blen);
    fprintf(stderr, "Initiated hash table.\n");

    while((l = kseq_read(seq)) >= 0) {
        tmp->bs_ptr = barcode_mem_view(seq);
        cp_view2buf(tmp->bs_ptr, tmp->key);
        HASH_FIND_STR(hash, tmp->key, tmp_hk);
        if(!tmp_hk) {
            tmp_hk = (HashKing_t *)malloc(sizeof(HashKing_t));
            tmp_hk->value = init_kfp(tmp->readlen);
            cp_view2buf(tmp->bs_ptr, tmp_hk->id);
            pushback_kseq(tmp_hk->value, seq, tmp->nuc_indices, tmp->blen);
            HASH_ADD_STR(hash, id, tmp_hk);
        }
        else {
            pushback_kseq(tmp_hk->value, seq, tmp->nuc_indices, tmp->blen);
        }
    }
    fprintf(stderr, "Loaded all fastq records into memory for meta-analysis. Now writing out to file!\n");
    HASH_ITER(hh, hash, current_entry, tmp_hk) {
        dmp_process_write(current_entry->value, out_handle, tmp->blen, tmp->buffers);
        destroy_kf(current_entry->value);
        free(current_entry->value);
        current_entry->value = NULL;
        HASH_DEL(hash, current_entry);
        free(current_entry);
    }
    fclose(out_handle);
    free(hash);
    kseq_destroy(seq);
    tmpvars_destroy(tmp);
}

