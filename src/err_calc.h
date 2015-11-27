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

void rate_calc(errcnt_t *e)
{
	uint64_t min_obs = min_obs;
	e->rates = (double ***)malloc(4 * sizeof(double **));
	e->qrates = (double **)malloc(4 * sizeof(double *));
	e->qdiffs = (int **)malloc(4 * sizeof(int *));
	uint64_t **qobs = (uint64_t **)malloc(4 * sizeof(uint64_t *));
	uint64_t **qerr = (uint64_t **)malloc(4 * sizeof(uint64_t *));
	double **qsum = (double **)malloc(4 * sizeof(double *));
	uint64_t **qcounts = (uint64_t **)malloc(4 * sizeof(uint64_t *));
	for(int i = 0; i < 4; ++i) {
		qcounts[i] = (uint64_t *)calloc(e->l, sizeof(uint64_t));
		qsum[i] = (double *)calloc(e->l, sizeof(double));
		qobs[i] = (uint64_t *)calloc(e->l,  sizeof(uint64_t));
		qerr[i] = (uint64_t *)calloc(e->l,  sizeof(uint64_t));
		e->qrates[i] = (double *)malloc(e->l * sizeof(double));
		e->rates[i] = (double **)malloc(nqscores * sizeof(double *));
		for(int j = 0; j < nqscores; ++j) {
			e->rates[i][j] = (double *)malloc(e->l * sizeof(double));
			for(uint32_t l = 0; l < e->l; ++l) {
				e->rates[i][j][l] = (e->obs[i][j][l] >= min_obs) ? (double)e->err[i][j][l] / e->obs[i][j][l] : DBL_MAX;
				qsum[i][l] += (-10. * log10(j + 2)) * e->obs[i][j][l];
				qcounts += e->obs[i][j][l];
				qobs[i][l] += e->obs[i][j][l];
				qerr[i][l] += e->err[i][j][l];
			}
		}
	}
	int i;
	uint64_t l;
	for(i = 0; i < 4; ++i) {
		for(l = 0; l < e->l; ++l)
			qsum[i][l] /= qcounts[i][l];
		free(qcounts[i]);
	}
	free(qcounts);

	for(i = 0; i < 4; ++i) {
		for(l = 0; l < e->l; ++l) {
			e->qrates[i][l] = (qobs[i][l]) ? (double)qerr[i][l] / qobs[i][l]: DBL_MAX;
			if(e->qrates[i][l] != DBL_MAX)
				e->qdiffs[i][l] = (int)(-10 * log10(e->qrates[i][l]) + 0.5) - (int)(-10 * log10(qsum[i][l]) + 0.5);
			else
				e->qdiffs[i][l] = 0; // Use ILMN-reported quality for these if not available
			// Difference between measured error rate and observed error rate
		}
		free(qobs[i]);
		free(qerr[i]);
	}
	e->final = (int ***)malloc(sizeof(int **) * 4);
	for(i = 0; i < 4; ++i) {
		e->final[i] = (int **)malloc(sizeof(int *) * nqscores);
		for(int j = 0; j < nqscores; ++j) {
			e->final[i][j] = (int *)malloc(sizeof(int) * e->l);
			for(l = 0; l < e->l; ++l)
				e->final[i][j][l] = (e->rates[i][j][l] == DBL_MAX) ?
						(int)(-10 * log10(e->rates[i][j][l]) + 0.5):
						j + 2 + e->qdiffs[i][l];
		}
	}
	free(qobs);
	free(qerr);
	// Now impute with the diffs for those which had insufficient observations.
	return;
}
// TODO: finish cleanup. Delete the qrates and qdiff arrays

errcnt_t *errcnt_init(size_t l) {
	errcnt_t *ret = (errcnt_t *)calloc(1, sizeof(errcnt_t));
	ret->obs = arr_init(l);
	ret->err = arr_init(l);
	ret->l = l;
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
