#include "include/kingfisher.h"
#include "dmp_interface.h"  // Contains the KingFisher_t type and KSEQ_INIT
#include "uthash_dmp_core.h"
#include "igamc_cephes.c"


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
    fprintf(stderr, "[omgz_core]: Now reading from file or handle %s.\n", strcmp(infname, "-") == 0 ? "stdin": infname);
    if(!outfname) out_handle = stdout;
    else {
        out_handle = fopen(outfname, "w");
    }
    gzFile fp = gzdopen(fileno(in_handle), "r");
    kseq_t *seq = kseq_init(fp);
    fprintf(stderr, "[omgz_core]: Opened file handles, initiated kseq parser.\n");
    // Initialized kseq
    int l = kseq_read(seq);
    if(l < 0) {
        fprintf(stderr, "[omgz_core]: Could not open fastq file (%s). Abort mission!\n",
                strcmp(infname, "-") == 0 ? "stdin": infname);
        exit(1);
    }
    char *bs_ptr = barcode_mem_view(seq);
    int blen = infer_barcode_length(bs_ptr);
    fprintf(stderr, "[omgz_core]: Barcode length: %i.\n", blen);
    tmpvars_t *tmp = init_tmpvars_p(bs_ptr, blen, seq->seq.l);
    // Start hash table
    HashKing_t *hash = NULL;
    HashKing_t *current_entry = (HashKing_t *)malloc(sizeof(HashKing_t));
    HashKing_t *tmp_hk = current_entry; // Save the pointer location for later comparison.
    cp_view2buf(bs_ptr, current_entry->id);
    current_entry->value = init_kfp(tmp->readlen);
    HASH_ADD_STR(hash, id, current_entry);
    pushback_kseq(current_entry->value, seq, tmp->nuc_indices, tmp->blen);

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
    fprintf(stderr, "[omgz_core]: Loaded all fastq records into memory for meta-analysis. Now writing out to file!\n");
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

