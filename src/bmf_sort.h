#ifndef BMF_SORT_H
#define BMF_SORT_H
#include <assert.h>
#include <errno.h>
#include <stdbool.h>
#include "htslib/khash.h"
#include "htslib/klist.h"
#include "htslib/ksort.h"
#include "htslib/kstring.h"
#include "dlib/cstr_util.h"
#include "dlib/io_util.h"
#include "dlib/misc_util.h"
#include "dlib/sort_util.h"
#include "dlib/bam_util.h"
#include "include/bam.h"
#include "include/sam_opts.h"


typedef bam1_t *bam1_p;
static inline int bam1_lt_ucs(const bam1_p a, const bam1_p b);
static inline int bam1_lt_bmf(const bam1_p a, const bam1_p b);

#ifdef __cplusplus
extern "C" {
#endif
    int sort_main(int argc, char **argv);
#ifdef __cplusplus
}
#endif


enum SORT_ORDER {
    QNAME,
    POS,
    BMF,
    UCS
};


#endif /* BMF_SORT_H */
