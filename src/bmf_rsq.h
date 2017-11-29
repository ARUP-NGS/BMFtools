#ifndef BMF_RSQ_H
#define BMF_RSQ_H
#include <assert.h>
#include "dlib/sort_util.h"
#include "dlib/bam_util.h"

#define STACK_START 128
#define READ_HD_LIMIT 6

namespace bmf {

enum cmpkey {
    POSITION,
    UNCLIPPED
};

CONST static inline int same_stack_pos_se(bam1_t *b, bam1_t *p)
{
    return (bmfsort_se_key(b) == bmfsort_se_key(p) &&
            b->core.l_qseq == p->core.l_qseq);
}

CONST static inline int same_stack_ucs_se(bam1_t *b, bam1_t *p)
{
    return (ucs_se_sort_key(b) == ucs_se_sort_key(p)  &&
             b->core.l_qseq == p->core.l_qseq);
}

struct StackFnPosSe {
    inline int operator()(bam1_t *b, bam1_t *p) const {return same_stack_pos_se(b, p);}
};

struct StackFnUcsSe {
    inline int operator()(bam1_t *b, bam1_t *p) const {return same_stack_ucs_se(b, p);}
};

CONST static inline int same_stack_pos(bam1_t *b, bam1_t *p)
{
    return (bmfsort_core_key(b) == bmfsort_core_key(p) &&
            bmfsort_mate_key(b) == bmfsort_mate_key(p));
}

CONST static inline int same_stack_ucs(bam1_t *b, bam1_t *p)
{
#if 0
    if(ucs_sort_core_key(b) == ucs_sort_core_key(p) &&
       ucs_sort_mate_key(b) == ucs_sort_mate_key(p) &&
       sort_rlen_key(b) == sort_rlen_key(p)) {
        assert(b->core.tid == p->core.tid);
        assert(b->core.mtid == p->core.mtid);
        assert(b->core.l_qseq == p->core.l_qseq);
        assert(bam_itag(b, "MU") == bam_itag(p, "MU"));
        assert(bam_itag(b, "LM") == bam_itag(p, "LM"));
        return 1;
    }
    return 0;
#else
    return (ucs_sort_core_key(b) == ucs_sort_core_key(p) &&
            ucs_sort_mate_key(b) == ucs_sort_mate_key(p));
#endif
}

struct StackFnPosPe {
    inline int operator()(bam1_t *b, bam1_t *p) const {return same_stack_pos(b, p);}
};

struct StackFnUcsPe {
    inline int operator()(bam1_t *b, bam1_t *p) const {return same_stack_ucs(b, p);}
};

CONST static inline int read_hd(const bam1_t *b, const bam1_t *p, const int lim=READ_HD_LIMIT)
{
    const uint8_t *const bseq = bam_get_seq(b);
    const uint8_t *const pseq = bam_get_seq(p);
    uint8_t bc, pc;
    int hd(0);
    for(int i(0); i < b->core.l_qseq; ++i) {
        bc = bam_seqi(bseq, i);
        pc = bam_seqi(pseq, i);
        if(bc != pc)
           if(bc != dlib::htseq::HTS_N)
               if(pc != dlib::htseq::HTS_N)
                   ++hd;
    }
    return hd;
}

CONST static inline int read_pass_hd(const bam1_t *b, const bam1_t *p, const int lim=READ_HD_LIMIT)
{
    const uint8_t *const bseq = bam_get_seq(b);
    const uint8_t *const pseq = bam_get_seq(p);
    uint8_t bc, pc;
    int hd(0);
    for(int i(0); i < b->core.l_qseq; ++i) {
        bc = bam_seqi(bseq, i);
        pc = bam_seqi(pseq, i);
        if(bc != pc &&
           bc != dlib::htseq::HTS_N &&
           pc != dlib::htseq::HTS_N &&
           ++hd > lim) // Note increment on hd.
            return 0;
    }
    return 1;
}

} /* namespace bmf */

#endif /* BMF_RSQ_H */
