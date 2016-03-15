#ifndef BMF_MARKRSQ_H
#define BMF_MARKRSQ_H
#include "dlib/bam_util.h"
#include "dlib/cstr_util.h"
#include "htslib/kstring.h"
#include <unistd.h>
#include <sys/stat.h>

namespace BMF {

    static inline void stack_bam1_cpy(bam1_t& dst, bam1_t *src) {
        uint8_t *data = static_cast<uint8_t *>(realloc(dst.data, src->m_data)); // save old data
        dst = *src; // Copy everything
        memcpy(data, src->data, src->l_data);
        dst.data = data;
    }

    struct bam1_stack {
        bam1_t *recs; // Pointer to array of recs. Access them by address.
        uint32_t n; // Number currently in use
        uint32_t m; // Number allocated
        void add_rec(bam1_t *b) {
            if(n + 1 > m) {
                const uint32_t old_m = m;
                kroundup32(m);
                recs = static_cast<bam1_t *>(realloc(recs, sizeof(bam1_t) * m));
                memset(recs + old_m, 0, (m - old_m) * sizeof(bam1_t)); // Set the new data members to null for realloc
            }
            stack_bam1_cpy(recs[n++], b);
        }
        void clear_stack_data() {
            // Since data is null or allocated, this is safe to do.
            while(m) free(recs[--m].data);
        }
        bam1_t *operator [](size_t index) {
            return recs + index;
        }
        bam1_stack(size_t start):
        recs(static_cast<bam1_t *>(calloc(start, sizeof(bam1_t) * start))),
        n(0),
        m(start)
        {
        }
        ~bam1_stack() {
            clear_stack_data();
            free(recs);
        }
    };
}
#endif
