#include "splitter.h"

#include <string.h>
#include "htslib/kstring.h"
#include "dlib/cstr_util.h"
#include "dlib/compiler_util.h"
#include "dlib/misc_util.h"
#include "lib/binner.h"

namespace BMF {

    void splitterhash_destroy(splitterhash_params_t *params)
    {
        /*
         * Is this freed by splitterhash_params?
        for(int i = 0; i < params->n; ++i) {
            LOG_DEBUG("i: %i.\n", i);
            cond_free(params->outfnames_r1[i]);
            cond_free(params->outfnames_r2[i]);
            cond_free(params->infnames_r1[i]);
            cond_free(params->infnames_r2[i]);
        }
        */
        cond_free(params->outfnames_r1);
        cond_free(params->outfnames_r2);
        cond_free(params->infnames_r1);
        cond_free(params->infnames_r2);
        cond_free(params);
    }

    void free_marksplit_settings(marksplit_settings_t settings)
    {
        cond_free(settings.tmp_basename);
        cond_free(settings.input_r1_path);
        cond_free(settings.input_r2_path);
        cond_free(settings.index_fq_path);
        cond_free(settings.rescaler);
        cond_free(settings.rescaler_path);
        cond_free(settings.homing_sequence);
        cond_free(settings.ffq_prefix);
    }

    splitterhash_params_t *init_splitterhash(marksplit_settings_t *settings_ptr, mark_splitter_t *splitter_ptr)
    {
        if(!settings_ptr) {
            LOG_EXIT("Settings pointer null. Abort!\n");
        }
        if(!settings_ptr->tmp_basename) {
            fprintf(stderr, "[E:%s] Output basename not set. Abort!\n", __func__);
            exit(EXIT_FAILURE);
        }
        if(!splitter_ptr) {
            fprintf(stderr, "[E:%s] Splitter pointer null. Abort!\n", __func__);
            exit(EXIT_FAILURE);
        }
        kstring_t ks = {0, 0, nullptr};
        splitterhash_params_t *ret = (splitterhash_params_t *)calloc(1, sizeof(splitterhash_params_t));
        ret->n = splitter_ptr->n_handles;
        if(settings_ptr->is_se) {
            ret->outfnames_r1 = (char **)malloc(ret->n * sizeof(char *));
            ret->infnames_r1 = (char **)malloc(ret->n * sizeof(char *));
            for(int i = 0; i < splitter_ptr->n_handles; ++i) {
                ret->infnames_r1[i] = splitter_ptr->fnames_r1[i];
                ks.l = 0;
                ksprintf(&ks, "%s.%i.dmp.fastq", settings_ptr->tmp_basename, i);
                ret->outfnames_r1[i] = ks_release(&ks);
            }
        } else {
            ret->outfnames_r1 = (char **)malloc(ret->n * sizeof(char *));
            ret->outfnames_r2 = (char **)malloc(ret->n * sizeof(char *));
            ret->infnames_r1 = (char **)malloc(ret->n * sizeof(char *));
            ret->infnames_r2 = (char **)malloc(ret->n * sizeof(char *));
            for(int i = 0; i < splitter_ptr->n_handles; ++i) {
                ret->infnames_r1[i] = splitter_ptr->fnames_r1[i];
                ret->infnames_r2[i] = splitter_ptr->fnames_r2[i]; // Does not allocate memory.  This is freed by mark_splitter_t!
                ks.l = 0;
                ksprintf(&ks, "%s.%i.R1.dmp.fastq", settings_ptr->tmp_basename, i);
                ret->outfnames_r1[i] = strdup(ks.s);
                ks.l = 0;
                ksprintf(&ks, "%s.%i.R2.dmp.fastq", settings_ptr->tmp_basename, i);
                ret->outfnames_r2[i] = strdup(ks.s);
            }
            free(ks.s);
        }
        return ret;
    }


    void splitter_destroy(mark_splitter_t *var)
    {
        for(int i = 0; i < var->n_handles; i++)
            cond_free(var->fnames_r1);
        if(var->fnames_r2)
            for(int i = 0; i < var->n_handles; ++i)
                cond_free(var->fnames_r2[i]);

        cond_free(var->tmp_out_handles_r1);
        cond_free(var->tmp_out_handles_r2);
        LOG_DEBUG("Freeing final var at %p.\n", (void *)var);
        cond_free(var);
    }


    mark_splitter_t init_splitter_pe(marksplit_settings_t* settings_ptr)
    {
        mark_splitter_t ret = {
            (gzFile *)malloc(settings_ptr->n_handles * sizeof(gzFile)), // tmp_out_handles_r1
            (gzFile *)malloc(settings_ptr->n_handles * sizeof(gzFile)), // tmp_out_handles_r2
            settings_ptr->n_nucs, // n_nucs
            (int)dlib::ipow(4, settings_ptr->n_nucs), // n_handles
            (char **)malloc(ret.n_handles * sizeof(char *)), // infnames_r1
            (char **)malloc(ret.n_handles * sizeof(char *))  // infnames_r2
        };
        kstring_t ks = {0, 0, nullptr};
        for (int i = 0; i < ret.n_handles; i++) {
            ks.l = 0;
            ksprintf(&ks, "%s.tmp.%i.R1.fastq", settings_ptr->tmp_basename, i);
            ret.fnames_r1[i] = dlib::kstrdup(&ks);
            ks.l = 0;
            ksprintf(&ks, "%s.tmp.%i.R2.fastq", settings_ptr->tmp_basename, i);
            ret.fnames_r2[i] = dlib::kstrdup(&ks);
            ret.tmp_out_handles_r1[i] = gzopen(ret.fnames_r1[i], settings_ptr->mode);
            ret.tmp_out_handles_r2[i] = gzopen(ret.fnames_r2[i], settings_ptr->mode);
        }
        return ret;
    }

    mark_splitter_t init_splitter_se(marksplit_settings_t* settings_ptr)
    {
        mark_splitter_t ret = {
            (gzFile *)calloc(settings_ptr->n_handles, sizeof(gzFile)), // tmp_out_handles_r1
            nullptr, // tmp_out_handles_r2
            settings_ptr->n_nucs, // n_nucs
            (int)dlib::ipow(4, settings_ptr->n_nucs), // n_handles
            (char **)calloc(ret.n_handles, sizeof(char *)), // infnames_r1
            nullptr  // infnames_r2
        };
        kstring_t ks = {0, 0, nullptr};
        for (int i = 0; i < ret.n_handles; i++) {
            ks.l = 0;
            ksprintf(&ks, "%s.tmp.%i.fastq", settings_ptr->tmp_basename, i);
            ret.fnames_r1[i] = dlib::kstrdup(&ks);
            ret.tmp_out_handles_r1[i] = gzopen(ret.fnames_r1[i], settings_ptr->mode);
        }
        free(ks.s);
        return ret;
    }

    mark_splitter_t init_splitter(marksplit_settings_t* settings_ptr)
    {
        if(settings_ptr->is_se) {
            return init_splitter_se(settings_ptr);
        } else {
            return init_splitter_pe(settings_ptr);
        }
    }

} /* namespace BMF */
