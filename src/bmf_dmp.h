#ifndef BMF_DMP_H
#define BMF_DMP_H

#include <assert.h>
#include <getopt.h>
#include <math.h>
#include <omp.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include "htslib/kstring.h"
#include "dlib/cstr_util.h"
#include "dlib/io_util.h"
#include "dlib/logging_util.h"
#include "dlib/math_util.h"
#include "dlib/mem_util.h"
#include "dlib/nix_util.h"
#include "lib/mseq.h"
#include "lib/binner.h"
#include "lib/kingfisher.h"
#include "bmf_hashdmp.h"

extern int64_t ipow(int32_t, int32_t);

typedef void (*hash_dmp_fn)(char *, char *, int);

#define RANDSTR_SIZE 20
#define DEFAULT_N_NUCS 4
#define DEFAULT_N_THREADS 4

namespace BMF {

    char test_hp_inline(char *barcode, int length, int threshold);
    void clean_homing_sequence(char *);
    void call_stdout(marksplit_settings_t *settings, splitterhash_params_t *params, char *ffq_r1, char *ffq_r2);
    void cat_fastqs(marksplit_settings_t *settings, splitterhash_params_t *params, char *ffq_r1, char *ffq_r2);
    void cat_fastqs_pe(marksplit_settings_t *settings, splitterhash_params_t *params, char *ffq_r1, char *ffq_r2);
    void cat_fastqs_se(marksplit_settings_t *settings, splitterhash_params_t *params, char *ffq_r1);
    void parallel_hash_dmp_core(marksplit_settings_t *settings, splitterhash_params_t *params, hash_dmp_fn func);
    void make_outfname(marksplit_settings_t *settings);
    void cleanup_hashdmp(marksplit_settings_t *settings, splitterhash_params_t *params);
    void check_rescaler(marksplit_settings_t *settings, int arr_size);

    /*
     * Returns 0 if a barcode is failed.
     * A barcode is failed by a homopolymer run >= threshold
     * or if an N is encountered.
     */
    CONST static inline int test_hp(char *barcode, int threshold)
    {
        assert(*barcode);
        int run = 0; char last = '\0';
        while(*barcode) {
            if(*barcode == 'N') return 0;
            if(*barcode == last) {
                if(++run == threshold) return 0;
            } else last = *barcode, run = 0;
            barcode++;
        }
        return 1;
    }


    /*
     * :param: settings [crms_settings_t, marksplit_settings_t] Settings struct in which to free the rescaler.
     * :return: void
     * This function supersedes free_rescaler_array by being type-generic.
     */

    #define cfree_rescaler(settings) \
        do {\
            if(settings.rescaler) {\
                int readlen##_settings = dlib::count_lines(settings.rescaler_path);\
                for(int i = 0; i < 2; ++i) {\
                    for(int j = 0; j < readlen##_settings; ++j) {\
                        for(int k = 0; k < nqscores; ++k) {\
                            cond_free(settings.rescaler[i][j][k]);\
                        }\
                        cond_free(settings.rescaler[i][j]);\
                    }\
                    cond_free(settings.rescaler[i]);\
                }\
                cond_free(settings.rescaler);\
            }\
        } while(0)


    static inline char *make_crms_outfname(char *fname)
    {
        return make_default_outfname(fname, ".crms.split");
    }


    static inline int nlen_homing_se(kseq_t *seq, marksplit_settings_t *settings_ptr, int default_len, int *pass_fail)
    {
        for(int i = settings_ptr->blen + settings_ptr->offset; i <= settings_ptr->max_blen; ++i) {
            if(memcmp(seq->seq.s + i, settings_ptr->homing_sequence, settings_ptr->homing_sequence_length) == 0) {
                *pass_fail = 1;
                return i + settings_ptr->homing_sequence_length;
            }
        }
        *pass_fail = 0;
        return default_len;
    }

    static inline int nlen_homing_default(kseq_t *seq1, kseq_t *seq2, marksplit_settings_t *settings_ptr, int default_len, int *pass_fail)
    {
        for(int i = settings_ptr->blen1_2 + settings_ptr->offset; i <= settings_ptr->max_blen; ++i) {
            if(!memcmp(seq1->seq.s + i, settings_ptr->homing_sequence, settings_ptr->homing_sequence_length)) {
                //LOG_DEBUG("Passed this one at %i.\n", i);
                *pass_fail = 1;
                return i + settings_ptr->homing_sequence_length;
            }
        }
        //LOG_INFO("Failed this on\n");
        *pass_fail = 0;
        return default_len;
    }

} /* namespace BMF */

#endif /* BMF_DMP_H */
