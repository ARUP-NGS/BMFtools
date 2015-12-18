#ifndef BMF_SORT_H
#define BMF_SORT_H
#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include <inttypes.h>
#include <getopt.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <zlib.h>
#include "bam.h"
#include "htslib/khash.h"
#include "htslib/klist.h"
#include "htslib/ksort.h"
#include "htslib/kstring.h"
#include "htslib/sam.h"
#include "flow_util.h"
#include "io_util.h"
#include "mem_util.h"
#include "sam.h"
#include "sam_opts.h"
#include "sort_util.h"
#include "cstr_util.h"

enum sort_order {
	SAMTOOLS,
	QNAME,
	BMF,
	UCS
};


#ifndef ABORT
#define ABORT(message) do {fprintf(stderr, message); exit(EXIT_FAILURE);} while(0)
#endif

static inline void add_unclipped_mate_starts(bam1_t *b1, bam1_t *b2) {
	uint32_t i, offset1 = 0, offset2 = 0, ucs1 = b1->core.pos, ucs2 = b2->core.pos;
	uint32_t *cigar1 = bam_get_cigar(b1);
	uint32_t *cigar2 = bam_get_cigar(b2);
	for(i = 0; i < b1->core.n_cigar; ++i) {
		if(!(cigar1[i]&0xf)) { // 'M' in cigar.
			break;
		}
		else {
			offset1 += cigar1[i] >> BAM_CIGAR_SHIFT;
		}
	}
	for(i = 0; i < b2->core.n_cigar; ++i) {
		if(!(cigar2[i]&0xf)) {
			break;
		}
		else {
			offset2 += cigar2[i] >> BAM_CIGAR_SHIFT;
		}
	}
	ucs1 += (b1->core.flag & BAM_FREVERSE) ? offset1: -1 * offset1;
	ucs2 += (b2->core.flag & BAM_FREVERSE) ? offset2: -1 * offset2;
	bam_aux_append(b2, "MU", 'I', sizeof(uint32_t), (uint8_t *)&ucs1);
	bam_aux_append(b1, "MU", 'I', sizeof(uint32_t), (uint8_t *)&ucs2);
	return;
}

#endif
