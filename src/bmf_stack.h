#ifndef BMF_STACK_H
#define BMF_STACK_H
#include "lib/UniqueObs.h"
#include "htslib/vcf.h"

int stack_main(int argc, char *argv[]);

#define DEFAULT_MAX_DEPTH (1 << 18)

struct stack_conf {
	float minFR; // Minimum fraction of family members agreed on base
	float minAF; // Minimum aligned fraction
	int max_depth;
	int minFM;
	uint32_t minFA;
	uint32_t minPV;
	int minMQ;
	int minCount;
	int minDuplex;
	int minOverlap;
	int skip_improper;
	uint32_t skip_flag; // Skip reads with any bits set to true
};

class stack_aux_t {
public:
	BamHandle *tumor;
	BamHandle *normal;
	vcfFile *ofp;
	bcf_hdr_t *vh;
	faidx_t *fai;
	khash_t(bed) *bed;
	int last_tid;
	char *ref_seq;
	stack_conf conf;
	stack_aux_t(char *tumor_path, char *normal_path, stack_conf _conf) {
		memset(this, 0, sizeof(*this));
		last_tid = -1;
		tumor = new BamHandle(tumor_path);
		normal = new BamHandle(normal_path);
		conf = _conf;
		if(!conf.max_depth) conf.max_depth = DEFAULT_MAX_DEPTH;
	}
	char get_ref_base(int tid, int pos) {
		int len;
		if(tid != last_tid) {
			if(ref_seq) free(ref_seq);
			ref_seq = fai_fetch(fai, tumor->header->target_name[tid], &len);
			last_tid = tid;
		}
		return ref_seq[pos];
	}
	~stack_aux_t() {
		if(ofp) vcf_close(ofp);
		if(bed) bed_destroy_hash((void *)bed);
		if(ref_seq) free(ref_seq);
		if(tumor) delete tumor;
		if(normal) delete normal;
	}
};


#endif /* BMF_STACK_H */
