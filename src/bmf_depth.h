#ifndef BMF_DEPTH_H
#define BMF_DEPTH_H
#include <zlib.h>
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "htslib/kstring.h"
#include "dlib/bam_util.h"
#include "dlib/cstr_util.h"
#include "dlib/mem_util.h"
#include "htslib/sam.h"
#include "sam_opts.h"

#include "htslib/kseq.h"

static inline int get_fm(bam1_t *b)
{
	uint8_t *data = bam_aux_get(b, "FM");
	return data ? bam_aux2i(data): 1;
}

#define DEFAULT_MAX_DEPTH (1 << 18)

int depth_main(int argc, char *argv[]);

#endif /* BMF_DEPTH_H */
