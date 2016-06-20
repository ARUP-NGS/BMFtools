#ifndef BMF_HASHDMP_H
#define BMF_HASHDMP_H
//#include "lib/mseq.h"
#include "lib/kingfisher.h"
extern "C" {
    #include "include/uthash.h"
}


#ifndef ifn_stream
#define ifn_stream(fname) ((fname) ? (fname): "stream")
#endif

namespace bmf {

    //KHASH_MAP_INIT_STR(dmp, kingfisher_t *)
    void hash_dmp_core(char *infname, char *outfname, int level);
    int hashdmp_main(int argc, char *argv[]);
    void stranded_hash_dmp_core(char *infname, char *outfname, int level);
    tmpvars_t *init_tmpvars_p(char *bs_ptr, int blen, int readlen);

    struct kingfisher_hash_t {
        UT_hash_handle hh;
        char id[MAX_BARCODE_LENGTH + 1];
        kingfisher_t *value;
    };


    CONST static inline int infer_barcode_length(char *bs_ptr)
    {
        char *const current = bs_ptr;
        for (;;) {
            switch(*bs_ptr++) {
            case '|': case '\0':
                return bs_ptr - current - 1;
            }
        }
        return -1; // This never happens.
    }

    static inline void cp_view2buf(char *view, char *buf)
    {
        for(;;) {
            switch(*view) {
                case '\0': case '|': *buf++ = '\0'; return;
                default: *buf++ = *view++;
            }
        }
    }

    static inline kingfisher_t *init_kfp(size_t readlen)
    {
        const size_t r5 = readlen * 5;
        kingfisher_t *ret = (kingfisher_t *)calloc(1, sizeof(kingfisher_t));
        ret->readlen = readlen;
        ret->max_phreds = (char *)malloc((r5) * sizeof(char));
        ret->nuc_counts = (uint16_t *)calloc(r5, sizeof(uint16_t));
        ret->phred_sums = (uint32_t *)calloc(r5, sizeof(uint32_t));
        memset(ret->max_phreds, '#', r5);
        ret->pass_fail = '1';
        return ret;
    }

    static inline void destroy_kf(kingfisher_t *kfp)
    {
        free(kfp->nuc_counts);
        free(kfp->phred_sums);
        free(kfp->max_phreds);
        free(kfp);
    }

    tmpvars_t *init_tmpvars_p(char *bs_ptr, int blen, int readlen);


    static inline void tmpvars_destroy(tmpvars_t *tmp)
    {
        free(tmp->buffers), free(tmp);
    }


    /* @func
     * Copies the barcode sequence from a fastq comment field into a buffer
     * :param: seq - [arg/kseq_t *] a filled-in kseq object.
     * :param: buf - a pre-allocated buffer or malloc'd char_ptr with enough space for the barcode and the null terminus.
     */
    #define cp_bs2buf(seq, buf) cp_view2buf(barcode_mem_view(seq), buf)


}

#endif /* BMF_HASHDMP_H */
