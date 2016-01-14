#ifndef BMF_ERR_H
#define BMF_ERR_H

#include <assert.h>
#include <float.h>
//#include "htslib/khash.h"
#include "htslib/faidx.h"
#include "dlib/bam_util.h"
#include "dlib/mem_util.h"
#include "dlib/bed_util.h"
#include "dlib/logging_util.h"
#include "lib/kingfisher.h"

#define DEFAULT_PADDING 50


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
	khash_t(bed) *bed; // parsed-in bed file hashmap. See dlib/bed_util.[ch] (http://github.com/NoSeatbelts/dlib).
	int minFM;
	int maxFM;
	int minMQ;
	int flag; // Filter flags. First use will simply be
} fullerr_t;

typedef struct obserr {
	uint64_t obs;
	uint64_t err;
} obserr_t;

KHASH_MAP_INIT_INT(obs, obserr_t)
KHASH_SET_INIT_INT(obs_union)

typedef struct fmerr {
	khash_t(obs) *hash1;
	khash_t(obs) *hash2;
	khash_t(bed) *bed;
    char *bedpath;
	char *refcontig;
	uint64_t flag;
	uint64_t nskipped;
	uint64_t nread;
    int minMQ;
} fmerr_t;

fmerr_t *fm_init(char *bedpath, bam_hdr_t *hdr, char *refcontig, int padding, int flag, int minMQ);
void fm_destroy(fmerr_t *fm);

enum err_flags {
	REQUIRE_DUPLEX = 1,
	REFUSE_DUPLEX = 2
};

fullerr_t *fullerr_init(size_t l, char *bedpath, bam_hdr_t *hdr,
        int padding, int minFM, int maxFM, int flag, int minMQ);
void fullerr_destroy(fullerr_t *e);

void err_core(char *fname, faidx_t *fai, fullerr_t *f, htsFormat *open_fmt);
void err_core_se(char *fname, faidx_t *fai, fullerr_t *f, htsFormat *open_fmt);


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
