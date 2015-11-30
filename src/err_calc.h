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
} fullerr_t;

#define nqscores 59uL
#define maxpv 3114

static inline int pv2ph(double pv)
{
	return (pv == 0.0) ? maxpv: (int)(-10. * log10(pv) + 0.5);
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

static const int bamseq2i[] = {-1, 0, 1, -1, 2, -1, -1, -1, 3};

void rate_calc(readerr_t *e);

readerr_t *readerr_init(size_t l) {
	readerr_t *ret = (readerr_t *)calloc(1, sizeof(readerr_t));
	arr3d_init(ret->obs, l, uint64_t);
	arr3d_init(ret->err, l, uint64_t);
	arr3d_init(ret->final, l, int);
	arr2d_init(ret->qdiffs, l, int);
	arr2d_init(ret->qpvsum, l, double);
	arr2d_init(ret->qobs, l, uint64_t);
	arr2d_init(ret->qerr, l, uint64_t);
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
	free(e);
	return;
}





#endif
