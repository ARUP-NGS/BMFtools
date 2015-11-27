#ifndef BED_UTIL_H
#define BED_UTIL_H

#include "khash.h"
#include "htslib/sam.h"
#include "htslib/vcf.h"
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <inttypes.h>

typedef struct region_set {
	uint64_t *intervals;
	uint64_t n;
} region_set_t;

KHASH_MAP_INIT_INT(bed, region_set_t)

int intcmp(const void *a, const void *b); // Compare intervals for sorting by start
void sort_bed(khash_t(bed) *bed);
khash_t(bed) *parse_bed_hash(char *path, bam_hdr_t *header, uint32_t padding);
void *bed_read(const char *fn);
void bed_destroy_hash(void *);
size_t get_nregions(khash_t(bed) *h);
static inline int bed_test(bam1_t *b, khash_t(bed) *h)
{
	khint_t k;
	if((k = kh_get(bed, h, b->core.tid)) == kh_end(h))
		return 0;
	const int32_t pos_plus_len = b->core.pos + b->core.l_qseq;
	for(uint64_t i = 0; i < kh_val(h, k).n; ++i) {
		if(pos_plus_len >= (int32_t)(kh_val(h, k).intervals[i] >> 32) && b->core.pos <= (int32_t)(kh_val(h, k).intervals[i]))
			return 1;
	}
	return 0;
}

static inline int vcf_bed_test(bcf1_t *b, khash_t(bed) *h)
{
	khint_t k;
	if((k = kh_get(bed, h, b->rid)) == kh_end(h))
		return 0;
	for(uint64_t i = 0; i < kh_val(h, k).n; ++i) {
		if(b->pos >= (int32_t)(kh_val(h, k).intervals[i] >> 32) && b->pos <= (int32_t)(kh_val(h, k).intervals[i]))
			return 1;
	}
	return 0;
}

#ifndef cond_free
#define cond_free(var) do {if(var) {free(var); var = NULL;}} while(0)
#endif

#endif /* BED_UTIL_H */
