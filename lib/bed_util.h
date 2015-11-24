#ifndef BED_UTIL_H
#define BED_UTIL_H

#include "khash.h"
#include "htslib/sam.h"
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>


typedef struct interval {
	uint32_t start;
	uint32_t end;
} interval_t;

typedef struct region_set {
	interval_t *intervals;
	uint64_t n;
} region_set_t;

KHASH_MAP_INIT_INT(bed, region_set_t)

int intcmp(const void *a, const void *b); // Compare intervals for sorting by start
void sort_bed(khash_t(bed) *bed);
khash_t(bed) *parse_bed(char *path, bam_hdr_t *header, int padding);
void *bed_read(char *fn);
void bed_destroy(void *);
static inline int bed_test(bam1_t *b, khash_t(bed) *h)
{
	khint_t k;
	k = kh_get(bed, h, b->core.tid);
	if(k == kh_end(h))
		return 0;
	uint32_t pos_plus_len = b->core.pos + b->core.l_qseq;
	for(uint64_t i = 0; i < kh_val(h, k).n; ++i) {
		if(pos_plus_len >= kh_val(h, k).intervals[i].start && b->core.pos <= kh_val(h, k).intervals[i].end)
			return 1;
	}
	return 0;
}

#ifndef cond_free
#define cond_free(var) do {if(var) {free(var); var = NULL;}} while(0)
#endif

#endif /* BED_UTIL_H */
