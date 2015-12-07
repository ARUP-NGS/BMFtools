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
#include "mem_util.h"
#include "seq_util.h"
#include "htslib/faidx.h"

typedef struct readerr {
	uint64_t ***obs;
	uint64_t ***err;
	uint64_t **qobs;
	uint64_t **qerr;
	double **qpvsum;
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
	char *refcontig;
} fullerr_t;


static inline int pv2ph(double pv)
{
	return (pv == 0.0) ? MAX_PV: (int)(-10. * log10(pv) + 0.5);
}

#define arr3d_init(var, l, type) \
	do {\
	var = (type ***)calloc(4, sizeof(type **));\
	for(int i_ = 0; i_ < 4; ++i_) {\
		var[i_] = (type **)calloc(nqscores, sizeof(type *));\
		for(int j_ = 0; j_ < nqscores; ++j_) {\
			var[i_][j_] = (type *)calloc(l, sizeof(type));\
		}\
	}} while(0)


#define arr2d_init(var, l, type) \
	do {\
	var = (type **)calloc(4, sizeof(type *));\
	for(int i_ = 0; i_ < 4; ++i_) {\
		var[i_] = (type *)calloc(l, sizeof(type));\
	}} while(0)

uint64_t ***a3d_u64(size_t l)
{
	uint64_t ***ret = (uint64_t ***)calloc(4, sizeof(uint64_t **));
	for(uint64_t i = 0; i < 4; ++i) {
		ret[i] = (uint64_t **)calloc(nqscores, sizeof(uint64_t *));
		for(uint64_t j = 0; j < nqscores; ++j)
			ret[i][j] = (uint64_t *)calloc(l, sizeof(uint64_t));
	}
	return ret;
}

double **a2d_double(size_t l)
{
	double **ret = (double **)calloc(4, sizeof(double *));
	for(int i = 0; i < 4; ++i)
		ret[i] = (double *)calloc(l, sizeof(double));
	return ret;
}

int **a2d_int(size_t l)
{
	int **ret = (int **)calloc(4, sizeof(int *));
	for(int i = 0; i < 4; ++i)
		ret[i] = (int *)calloc(l, sizeof(int));
	return ret;
}

uint64_t **a2d_u64(size_t l)
{
	uint64_t **ret = (uint64_t **)calloc(4, sizeof(uint64_t *));
	for(uint64_t i = 0; i < 4; ++i)
		ret[i] = (uint64_t *)calloc(l, sizeof(uint64_t));
	return ret;
}


int ***a3d_int(size_t l)
{
	int ***ret = (int ***)calloc(4, sizeof(int **));
	for(int i = 0; i < 4; ++i) {
		ret[i] = (int **)calloc(nqscores, sizeof(int *));
		for(int j = 0; j < nqscores; ++j) {
			ret[i][j] = (int *)calloc(l, sizeof(int));
		}
	}
	return ret;
}

static const int bamseq2i[] = {-1, 0, 1, -1, 2, -1, -1, -1, 3};

void rate_calc(readerr_t *e);

readerr_t *readerr_init(size_t l) {
	readerr_t *ret = (readerr_t *)calloc(1, sizeof(readerr_t));
	ret->obs = a3d_u64(l);
	ret->err = a3d_u64(l);
	ret->final = a3d_int(l);
	ret->qdiffs = a2d_int(l);
	ret->qpvsum = a2d_double(l);
	ret->qobs = a2d_u64(l);
	ret->qerr = a2d_u64(l);
	ret->l = l;
	return ret;
}

fullerr_t *fullerr_init(size_t l) {
	fullerr_t *ret = (fullerr_t *)calloc(1, sizeof(fullerr_t));
	ret->l = l;
	ret->r1 = readerr_init(l);
	ret->r2 = readerr_init(l);
	return ret;
}

void readerr_destroy(readerr_t *e);
void fullerr_destroy(fullerr_t *e) {
#if !NDEBUG
	fprintf(stderr, "Beginning fullerr_destroy\n");
#endif
	if(e->r1)
		readerr_destroy(e->r1), e->r1 = NULL;
	if(e->r2)
		readerr_destroy(e->r2), e->r2 = NULL;
	if(e->refcontig) free(e->refcontig), e->refcontig = NULL;
	free(e);
	return;
}





#endif
