#ifndef BMF_SORT_H
#define BMF_SORT_H
#include <assert.h>
#include <errno.h>
#include <inttypes.h>
#include <stdbool.h>
#include "htslib/khash.h"
#include "htslib/klist.h"
#include "htslib/ksort.h"
#include "htslib/kstring.h"
#include "dlib/cstr_util.h"
#include "dlib/io_util.h"
#include "dlib/mem_util.h"
#include "dlib/sort_util.h"
#include "dlib/bam_util.h"
#include "include/bam.h"
#include "include/sam_opts.h"

#ifdef __cplusplus
extern "C" {
#endif
    int sort_main(int argc, char **argv);
    extern char *rand_string(char *, size_t size);
#ifdef __cplusplus
}
#endif

enum sort_order {
    SAMTOOLS,
    QNAME,
    BMF_POS,
    UCS
};


#endif /* BMF_SORT_H */
