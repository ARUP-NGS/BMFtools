#ifndef BMF_STACK_H
#define BMF_STACK_H
#include "lib/UniqueObs.h"
#include "htslib/vcf.h"

int stack_main(int argc, char *argv[]);

#define DEFAULT_MAX_DEPTH (1 << 18)

class aux_t {
public:
	samFile *fp;
	hts_itr_t *iter;
	hts_idx_t *bam_index;
	bam_hdr_t *bh;
	vcfFile *ofp;
	bcf_hdr_t *vh;
	khash_t(bed) *bed;
	float minFR; // Minimum fraction of family members agreed on base
	float minAF; // Minimum aligned fraction
	int max_depth;
	int minFM;
	uint32_t minFA;
	uint32_t minPV;
	uint32_t minMQ;
	int minCount;
	int minDuplex;
	int minOverlap;
	int skip_improper;
	uint32_t skip_flag; // Skip reads with any bits set to true
	aux_t() {
		memset(this, 0, sizeof(*this));
		max_depth = DEFAULT_MAX_DEPTH;
	}
	~aux_t() {
		if(fp) sam_close(fp);
		if(iter) hts_itr_destroy(iter);
		if(bam_index) hts_idx_destroy(bam_index);
		if(bh) bam_hdr_destroy(bh);
		if(ofp) vcf_close(ofp);
		if(bed) bed_destroy_hash((void *)bed);
	}
};


#endif /* BMF_STACK_H */
