#ifndef ERR_CALC_H
#define ERR_CALC_H

#include <stdint.h>
#include <stddef.h>
#include <inttypes.h>
#include <assert.h>
#include <float.h>
#include <math.h>
#include "khash.h"
#include "pair_util.h"
#include "htslib/faidx.h"

typedef struct readerr {
	uint64_t ***obs;
	uint64_t ***err;
	uint64_t ***qcounts;
	double ***rates;
	double **qrates;
	int **qdiffs;
	int ***final;
	size_t l; // Read length
} readerr_t;

typedef struct fullerr {
	uint64_t nread; // Number of records read
	uint64_t nskipped; // Number of records read
	readerr_t *r1;
	readerr_t *r2;
	size_t l;
} fullerr_t;

#define nqscores 49uL

static int bamseq2i[8] = {-1, 0, 1, -1, 2, -1, -1, 3};
static uint64_t min_obs = 1000;

static inline int pv2ph(double pv)
{
	return (int)(-10. * log10(pv) + 0.5);
}

uint64_t ***arr_init(size_t l) {
	uint64_t ***ret = (uint64_t ***)calloc(4, sizeof(uint64_t **));
	for(int i = 0; i < 4; ++i) {
		uint64_t **q_arrs = (uint64_t **)calloc(nqscores, sizeof(uint64_t *));
		for(int j = 0; j < nqscores; ++j) {
			uint64_t *arr_rates = (uint64_t *)calloc(l, sizeof(uint64_t));
			q_arrs[j] = arr_rates;
		}
		ret[i] = q_arrs;
	}
	return ret;
}

void rate_calc(readerr_t *e);

readerr_t *readerr_init(size_t l) {
	readerr_t *ret = (readerr_t *)calloc(1, sizeof(readerr_t));
	ret->obs = arr_init(l);
	ret->err = arr_init(l);
	ret->l = l;
	return ret;
}

fullerr_t *fullerr_init(size_t l) {
	fullerr_t *ret = (fullerr_t *)calloc(1, sizeof(fullerr_t));
	ret->r1 = readerr_init(l);
	ret->r2 = readerr_init(l);
	ret->l = l;
	return ret;
}
void readerr_destroy(readerr_t *e);
void fullerr_destroy(fullerr_t *e) {
	if(e->r1)
		readerr_destroy(e->r1), e->r1 = NULL;
	if(e->r2)
		readerr_destroy(e->r2), e->r2 = NULL;
	free(e);
	return;
}





#endif
