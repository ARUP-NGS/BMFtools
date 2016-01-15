#ifndef BMF_FAMSTATS_H
#define BMF_FAMSTATS_H

#include <getopt.h>

#include "htslib/sam.h"
#include "htslib/ksort.h"
#include "htslib/khash.h"
#include "dlib/bam_util.h"
#include "dlib/bed_util.h"
#include "dlib/compiler_util.h"


KHASH_MAP_INIT_INT64(fm, uint64_t)
KHASH_MAP_INIT_INT64(rc, uint64_t)


#ifndef cond_free
#define cond_free(var) do {if(var) {free(var); var = NULL;}} while(0)
#endif

typedef struct famstats {
	uint64_t n_pass;
	uint64_t n_fail;
	uint64_t allfm_sum;
	uint64_t allfm_counts;
	uint64_t allrc_sum;
	uint64_t realfm_sum;
	uint64_t realfm_counts;
	uint64_t realrc_sum;
	uint64_t dr_sum;
	uint64_t dr_counts;
	uint64_t dr_rc_sum;
	double dr_rc_frac_sum;
	khash_t(fm) *fm;
	khash_t(rc) *rc;
	khiter_t ki;
	int khr;
	uint8_t *data;
} famstats_t;

typedef struct famstat_settings {
	uint64_t notification_interval;
	uint32_t minMQ;
	uint32_t minFM;
} famstat_settings_t;
#endif /* BMF_FAMSTATS_H */
