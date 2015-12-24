#ifndef BMF_ERR_H
#define BMF_ERR_H

#include <assert.h>
#include <float.h>
#include "htslib/khash.h"
#include "htslib/faidx.h"
#include "dlib/bam_util.h"
#include "dlib/mem_util.h"
#include "lib/kingfisher.h"


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

readerr_t *readerr_init(size_t l);
void readerr_destroy(readerr_t *e);
void rate_calc(readerr_t *e);


typedef struct fullerr {
	uint64_t nread; // Number of records read
	uint64_t nskipped; // Number of records read
	readerr_t *r1;
	readerr_t *r2;
	size_t l;
	char *refcontig;
} fullerr_t;

fullerr_t *fullerr_init(size_t l);
void fullerr_destroy(fullerr_t *e);


static inline int pv2ph(double pv)
{
	return (pv > 0.0) ? (int)(-10. * log10(pv) + 0.5): MAX_PV;
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


#endif /* BMF_ERR_H */
