#include "include/kseq_dec.h"


typedef struct tmpvars {
    char *bs_ptr;
    int blen;
    int readlen;
    int nuc_indices[2];
    char key[MAX_BARCODE_LENGTH + 1];
    int l; // For holding ret value for seq.
    tmpbuffers_t *buffers;
} tmpvars_t;


typedef struct HashKing {
    UT_hash_handle hh;
    char id[MAX_BARCODE_LENGTH + 1];
    KingFisher_t *value;
} HashKing_t;

typedef struct splitterhash_params {
    char **infnames_r1;
    char **infnames_r2;
    char **outfnames_r1;
    char **outfnames_r2;
    int n; // Number of infnames and outfnames
    int paired; // 1 if paired, 0 if single-end
} splitterhash_params_t;


/*
 * :param: seq - [arg/kseq_t *] a filled-in kseq object.
 * :param: buf - a pre-allocated buffer or malloc'd char_ptr with enough space for the barcode and the null terminus.
 * :returns:
 */
inline void cp_bs2buf(kseq_t *seq, char *buf)
{
    char *view = barcode_mem_view(seq);
    int blen = 0;
    while(view[blen] != '\0' && view[blen] != '|') {
        ++blen;
    }
    memcpy(buf, view, blen);
    buf[blen] = '\0';
    return;
}


inline splitterhash_params_t *init_splitterhash(mssi_settings_t *settings_ptr, mark_splitter_t *splitter_ptr)
{
    char tmp_buffer [METASYNTACTIC_FNAME_BUFLEN];
    splitterhash_params_t *ret = (splitterhash_params_t *)malloc(sizeof(splitterhash_params_t));
    ret->n = splitter_ptr->n_handles;
    ret->outfnames_r1 = (char **)malloc(ret->n * sizeof(char *));
    for(int i = 0; i < splitter_ptr->n_handles; ++i) {
        ret->infnames_r1[i] = splitter_ptr->fnames_r1[i];
        ret->infnames_r2[i] = splitter_ptr->fnames_r2[i]; // Does not allocate memory.  This is freed by mark_splitter_t!
        sprintf(tmp_buffer, "%s.%i.R1.dmp.fastq", settings_ptr->output_basename, i);
        ret->outfnames_r1[i] = strdup(tmp_buffer);
        sprintf(tmp_buffer, "%s.%i.R2.dmp.fastq", settings_ptr->output_basename, i);
        ret->outfnames_r2[i] = strdup(tmp_buffer);
    }
}

inline void splitterhash_destroy(splitterhash_params_t *params)
{
    for(int i = 0; i < params->n; ++i) {
        free(params->outfnames_r1[i]);
        free(params->outfnames_r2[i]);
    }
    free(params->outfnames_r1);
    free(params->outfnames_r2);
    free(params);
    params = NULL;
    return;
}
