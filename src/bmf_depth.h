#ifndef BMF_DEPTH_H
#define BMF_DEPTH_H
#include <assert.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <zlib.h>
#include <vector>
#include <algorithm>
#include "htslib/kstring.h"
#include "dlib/bed_util.h"
#include "dlib/bam_util.h"
#include "dlib/cstr_util.h"
#include "dlib/mem_util.h"
#include "htslib/sam.h"
#include "sam_opts.h"

#include "htslib/kseq.h"

#define DEFAULT_MAX_DEPTH (1 << 18)

KHASH_MAP_INIT_INT(depth, uint64_t)

int depth_main(int argc, char *argv[]);

#endif /* BMF_DEPTH_H */
