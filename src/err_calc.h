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

typedef struct errcnt {
	uint64_t ***obs;
	uint64_t ***err;
	uint64_t ***qcounts;
	double ***rates;
	double **qrates;
	int **qdiffs;
	int ***final;
	size_t l; // Read length
	uint64_t nread; // Number of records read
	uint64_t nskipped; // Number of records read
} errcnt_t;

#define nqscores 43

static int bamseq2i[8] = {-1, 0, 1, -1, 2, -1, -1, 3};
static const uint64_t min_obs = 1000;

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

void rate_calc(errcnt_t *e);

errcnt_t *errcnt_init(size_t l) {
	errcnt_t *ret = (errcnt_t *)calloc(1, sizeof(errcnt_t));
	ret->obs = arr_init(l);
	ret->err = arr_init(l);
	ret->l = l;
	return ret;
}


void errcnt_destroy(errcnt_t *e);


#endif
