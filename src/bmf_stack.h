#ifndef BMF_STACK_H
#define BMF_STACK_H
#include "lib/UniqueObs.h"
#include "htslib/vcf.h"

int stack_main(int argc, char *argv[]);

#define DEFAULT_MAX_DEPTH (1 << 18)

class stack_aux_t {
public:
	samFile *fp;
	hts_itr_t *iter;
	hts_idx_t *bam_index;
	bam_hdr_t *bh;
	vcfFile *ofp;
	bcf_hdr_t *vh;
	faidx_t *fai;
	khash_t(bed) *bed;
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
	int last_tid;
	char *ref_seq;
	uint32_t skip_flag; // Skip reads with any bits set to true
	stack_aux_t() {
		memset(this, 0, sizeof(*this));
		max_depth = DEFAULT_MAX_DEPTH;
		last_tid = -1;
	}
	char get_ref_base(int tid, int pos) {
		int len;
		if(tid != last_tid) {
			if(ref_seq) free(ref_seq);
			ref_seq = fai_fetch(fai, bh->target_name[tid], &len);
			last_tid = tid;
		}
		return ref_seq[pos];
	}
	~stack_aux_t() {
		if(fp) sam_close(fp);
		if(iter) hts_itr_destroy(iter);
		if(bam_index) hts_idx_destroy(bam_index);
		if(bh) bam_hdr_destroy(bh);
		if(ofp) vcf_close(ofp);
		if(bed) bed_destroy_hash((void *)bed);
		if(ref_seq) free(ref_seq);
	}
};


#endif /* BMF_STACK_H */
