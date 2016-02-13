#ifndef BMF_ERR_H
#define BMF_ERR_H

#include <assert.h>
#include <float.h>
#include "htslib/faidx.h"
#include "dlib/bam_util.h"
#include "dlib/mem_util.h"
#include "dlib/bed_util.h"
#include "dlib/logging_util.h"
#include "lib/kingfisher.h"
#include <string>
#include <vector>


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
	uint32_t minPV;
	int flag; // Filter flags. First use will simply be
} fullerr_t;

typedef struct obserr {
	uint64_t obs;
	uint64_t err;
} obserr_t;

class RegionErr {
	obserr_t counts;
public:
	std::string name;
	inline void inc_err() {
		++counts.err;
	}
	inline void inc_obs() {
		++counts.obs;
	}
	RegionErr(region_set_t set, int i);
	inline void write_report(FILE *fp) {
		fprintf(fp, "%s\t", name.c_str());
		if(counts.obs) {
			fprintf(fp, "%s\t%f\t%lu\t%lu\n", name.c_str(), ((double)counts.err * 100.) / counts.obs,
					counts.err, counts.obs);
		} else {
			fprintf(fp, "%s\t-nan\t0\t0\n", name.c_str());
		}
	}
};

typedef struct cycle_err {
	int32_t minMQ;
	int32_t rlen;
	obserr_t *r1;
	obserr_t *r2;
	uint64_t nskipped;
	uint64_t nread;
	char *refcontig;
	char *bedpath;
	int flag;
	khash_t(bed) *bed; // parsed-in bed file hashmap. See dlib/bed_util.[ch] (http://github.com/NoSeatbelts/dlib).
} cycle_err_t;

class RegionExpedition {
	uint32_t padding;
	std::string bedpath;
public:
	khash_t(bed) *bed;
	faidx_t *fai; // Does not own fai!
	std::vector<RegionErr> region_counts;
	samFile *fp;
	bam_hdr_t *hdr;
	hts_itr_t *iter;
	hts_idx_t *bam_index;
	int32_t minMQ;
	int32_t minFM;
	int32_t requireFP;
	int32_t max_depth;
	RegionExpedition(char *bampath, char *_bedpath, faidx_t *_fai, int32_t _minMQ=0, uint32_t _padding=DEFAULT_PADDING,
			int32_t _minFM=0, int32_t _requireFP=0, int _max_depth=262144) {
		fp = sam_open(bampath, "r");
		hdr = sam_hdr_read(fp);
		if(!fp || !hdr) LOG_ERROR("Could not open input sam file %s. Abort!\n", bampath);
		padding = _padding;
		bed = parse_bed_hash(_bedpath, hdr, _padding);
		bedpath = std::string(_bedpath);
		minMQ = _minMQ;
		minFM = _minFM;
		max_depth = _max_depth;
		requireFP = _requireFP;
		fai = _fai;
		region_counts = std::vector<RegionErr>();
		iter = NULL;
		bam_index= sam_index_load(fp, fp->fn);
		if(!bam_index) LOG_ERROR("Could not read bam index for sam file %s. Abort!\n", fp->fn);
	}
	~RegionExpedition() {
		if(bed) bed_destroy_hash((void *)bed);
		if(iter) hts_itr_destroy(iter);
		if(bam_index) hts_idx_destroy(bam_index);
		if(hdr) bam_hdr_destroy(hdr);
		if(fp) sam_close(fp);
	}
	char *get_bampath() {
		return fp ? fp->fn: NULL;
	}
};

cycle_err_t *cycle_init(char *bedpath, bam_hdr_t *hdr, char *refcontig, int padding, int minMQ, int rlen, int flag);
void cycle_destroy(cycle_err_t *c);
extern void check_bam_tag_exit(char *bampath, const char *tag);

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
	unsigned minMQ;
	unsigned minPV;
} fmerr_t;

fmerr_t *fm_init(char *bedpath, bam_hdr_t *hdr, char *refcontig, int padding, int flag, unsigned minMQ, uint32_t minPV);
void fm_destroy(fmerr_t *fm);

enum err_flags {
	REQUIRE_DUPLEX = 1,
	REFUSE_DUPLEX = 2,
	REQUIRE_PROPER = 4,
	REQUIRE_FP_PASS = 8
};

fullerr_t *fullerr_init(size_t l, char *bedpath, bam_hdr_t *hdr,
		int padding, int minFM, int maxFM, int flag, int minMQ, uint32_t minPV);
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
	for(unsigned i_ = 0; i_ < 4u; ++i_) {\
		var[i_] = (type **)calloc(nqscores, sizeof(type *));\
		for(unsigned j_ = 0; j_ < nqscores; ++j_) {\
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


#define cycle_loop(ce, b, seq, cigar, last_tid, ref, hdr, pos, is_rev, length, i, rc, fc, ind, s, cycle, obsarr, rlen, fpdata)\
	do {\
		LOG_DEBUG("Got into cycle loop.\n");\
		if(++ce->nread % 1000000 == 0) {\
			LOG_INFO("Records read: %lu.\n", ce->nread);\
		}\
		if((b->core.flag & 772) || \
			b->core.qual < ce->minMQ || \
			(ce->refcontig && tid_to_study != b->core.tid) ||\
			((ce->flag & REQUIRE_FP_PASS) && ((fpdata = bam_aux_get(b, "FP")) != NULL) && bam_aux2i(fpdata) == 0) || \
			((ce->flag & REQUIRE_PROPER) && (!(b->core.flag & BAM_FPROPER_PAIR))) ||\
			(ce->bed && bed_test(b, ce->bed) == 0) /* Outside of region */) {\
			++ce->nskipped;\
			LOG_DEBUG("Skipped record with name %s.\n", bam_get_qname(b));\
			continue;\
		}\
		seq = (uint8_t *)bam_get_seq(b);\
		cigar = bam_get_cigar(b);\
		if(b->core.tid != last_tid) {\
			cond_free(ref);\
			LOG_DEBUG("Loading ref sequence for contig with name %s.\n", hdr->target_name[b->core.tid]);\
			ref = fai_fetch(fai, hdr->target_name[b->core.tid], &reflen);\
			if(!ref) {\
				LOG_ERROR("Failed to load ref sequence for contig '%s'. Abort!\n", hdr->target_name[b->core.tid]);\
			}\
			last_tid = b->core.tid;\
		}\
		pos = b->core.pos;\
		is_rev = b->core.flag & BAM_FREVERSE;\
		obsarr = (b->core.flag & BAM_FREAD1) ? ce->r1: ce->r2;\
		LOG_DEBUG("obsarr pointer: %p.\n", (void *)obsarr);\
		for(i = 0, rc = 0, fc = 0; i < b->core.n_cigar; ++i) {\
			length = bam_cigar_oplen(cigar[i]);\
			switch(bam_cigar_op(cigar[i])) {\
			case BAM_CMATCH:\
			case BAM_CEQUAL:\
			case BAM_CDIFF:\
				for(ind = 0; ind < length; ++ind) {\
					s = bam_seqi(seq, ind + rc);\
					if(s == HTS_N || ref[pos + fc + ind] == 'N') continue;\
					cycle = is_rev ? rlen - 1 - ind - rc: ind + rc;\
					++obsarr[cycle].obs;\
					if(seq_nt16_table[(int8_t)ref[pos + fc + ind]] != s) ++obsarr[cycle].err;\
				}\
				rc += length; fc += length;\
				break;\
			case BAM_CSOFT_CLIP:\
			case BAM_CHARD_CLIP:\
			case BAM_CINS:\
				rc += length;\
				break;\
			case BAM_CREF_SKIP:\
			case BAM_CDEL:\
				fc += length;\
				break;\
			}\
		}\
	} while(0)

#endif /* BMF_ERR_H */
