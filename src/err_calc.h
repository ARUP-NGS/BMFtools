#ifndef ERR_CALC_H
#define ERR_CALC_H

#include <stdint.h>
#include <stddef.h>
#include <inttypes.h>
#include <assert.h>
#include <float.h>
#include "khash.h"
#include "pair_util.h"
#include "htslib/faidx.h"

typedef struct errcnt {
	uint64_t ***obs;
	uint64_t ***err;
	double ***rates;
	size_t l; // Read length
	uint64_t nread; // Number of records read
	uint64_t nskipped; // Number of records read
} errcnt_t;

#define nqscores 59

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

void rate_init(errcnt_t *e)
{
	e->rates = (double ***)malloc(4 * sizeof(double **));
	for(int i = 0; i < 4; ++i) {
		e->rates[i] = (double **)malloc(nqscores * sizeof(double *));
		for(int j = 0; j < nqscores; ++j) {
			e->rates[i][j] = (double *)malloc(e->l * sizeof(double));
			for(uint32_t l = 0; l < e->l; ++l) {
				e->rates[i][j][l] = (double)e->err[i][j][l] / e->obs[i][j][l];
			}
		}
	}
	return;
}

errcnt_t *errcnt_init(size_t l) {
	errcnt_t *ret = (errcnt_t *)calloc(1, sizeof(errcnt_t));
	ret->obs = arr_init(l);
	ret->err = arr_init(l);
	ret->l = l;
	rate_init(ret);
	return ret;
}


void errcnt_destroy(errcnt_t *e) {
	for(int i = 0; i < 4; ++i) {
		for(int j = 0; j < nqscores; ++j) {
			free(e->obs[i][j]);
			free(e->err[i][j]);
		}
		free(e->obs[i]);
		free(e->err[i]);
	}
	if(e->rates) {
		for(int i = 0; i < 4; ++i) {
			for(int j = 0; j < nqscores; ++j) {
				free(e->rates[i][j]);
			}
			free(e->rates[i]);
		}
		free(e->rates);
		e->rates = NULL;
	}
	free(e->obs);
	free(e->err);
	free(e);
}

#endif
