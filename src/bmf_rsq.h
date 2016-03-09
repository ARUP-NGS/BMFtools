#ifndef BMF_RSQ_H
#define BMF_RSQ_H
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <zlib.h>
#include <tgmath.h>
#include <unordered_map>
#include <algorithm>
#include <functional>
#include "htslib/sam.h"
#include "include/sam_opts.h"
#include "include/bam.h" // for bam_get_library
#include "include/igamc_cephes.h" /// for igamc
#include "dlib/cstr_util.h"
#include "dlib/sort_util.h"
#include "dlib/bam_util.h"

#define STACK_START 128
#define READ_HD_LIMIT 6

namespace BMF {

    struct tmp_stack_t {
        size_t n, max;
        bam1_t **a;
    };

    static inline void stack_insert(tmp_stack_t *stack, bam1_t *b)
    {
        if (stack->n == stack->max) {
            stack->max = stack->max? stack->max<<1 : 0x10000;
            stack->a = (bam1_t**)realloc(stack->a, sizeof(bam1_t*) * stack->max);
        }
        stack->a[stack->n++] = bam_dup1(b);
    }


    enum cmpkey {
        POSITION,
        UNCLIPPED
    };

    CONST static inline int same_stack_pos_se(bam1_t *b, bam1_t *p)
    {
        return bmfsort_se_key(b) == bmfsort_se_key(p);
    }

    CONST static inline int same_stack_ucs_se(bam1_t *b, bam1_t *p)
    {
        return ucs_se_sort_key(b) == ucs_se_sort_key(p);
    }

    CONST static inline int same_stack_pos(bam1_t *b, bam1_t *p)
    {
        return (bmfsort_core_key(b) == bmfsort_core_key(p) &&
                bmfsort_mate_key(b) == bmfsort_mate_key(p));
    }

    CONST static inline int same_stack_ucs(bam1_t *b, bam1_t *p)
    {
    #if !NDEBUG
        if(ucs_sort_core_key(b) == ucs_sort_core_key(p) &&
                ucs_sort_mate_key(b) == ucs_sort_mate_key(p)) {
            assert(b->core.tid == p->core.tid);
            assert(b->core.mtid == p->core.mtid);
            assert(b->core.mtid == p->core.mtid);
            assert(bam_itag(b, "MU") == bam_itag(p, "MU"));
            assert(bam_itag(b, "SU") == bam_itag(p, "SU"));
            return 1;
        }
        return 0;
    #else
        return (ucs_sort_core_key(b) == ucs_sort_core_key(p) &&
                ucs_sort_mate_key(b) == ucs_sort_mate_key(p));
    #endif
    }


    CONST static inline int read_pass_hd(bam1_t *b, bam1_t *p, const int lim=READ_HD_LIMIT)
    {
        const uint8_t *const bseq = bam_get_seq(b);
        const uint8_t *const pseq = bam_get_seq(p);
        uint8_t bc, pc;
        int hd = 0;
        for(int i = 0; i < b->core.l_qseq; ++i) {
            bc = bam_seqi(bseq, i);
            pc = bam_seqi(pseq, i);
            if(bc != pc && bc != dlib::htseq::HTS_N && pc != dlib::htseq::HTS_N && ++hd > lim)
                return 0;
        }
        return 1;
    }

} /* namespace BMF */

#endif /* BMF_RSQ_H */
