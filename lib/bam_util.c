#include "bam_util.h"

bam_plp_t bam_plp_maxcnt_init(bam_plp_auto_f func, void *data, int maxcnt)
{
	bam_plp_t iter = bam_plp_init(func, data);
	bam_plp_set_maxcnt(iter, maxcnt);
    return iter;
}

void abstract_single_data(samFile *in, bam_hdr_t *hdr, samFile *out, single_aux function, void *data)
{
	bam1_t *b = bam_init1();
	while (LIKELY(sam_read1(in, hdr, b) >= 0))
		function(b, data), sam_write1(out, hdr, b);
	bam_destroy1(b);
}

void abstract_single_iter(samFile *in, bam_hdr_t *hdr, samFile *out, single_fn function)
{
	bam1_t *b = bam_init1();
	while (LIKELY(sam_read1(in, hdr, b) >= 0))
		function(b), sam_write1(out, hdr, b);
	bam_destroy1(b);
}

void abstract_single_filter(samFile *in, bam_hdr_t *hdr, samFile *out, single_aux_check function, void *data)
{
	bam1_t *b;
	b = bam_init1();
	while (LIKELY(sam_read1(in, hdr, b) >= 0)) {
		if(function(b, data))
			continue;
		sam_write1(out, hdr, b);
	}
	bam_destroy1(b);
}

void abstract_pair_iter(samFile *in, bam_hdr_t *hdr, samFile *ofp, pair_fn function)
{
	bam1_t *b = bam_init1(), *b1 = bam_init1();
	while (LIKELY(sam_read1(in, hdr, b) >= 0)) {
		bam1_core_t *c = &b->core;
		if(c->flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY))
			continue;
		if(c->flag & BAM_FREAD1) {
			b1 = bam_copy1(b1, b);
			continue; // b w
		}
		function(b1, b);
		sam_write1(ofp, hdr, b), sam_write1(ofp, hdr, b1);
	}
	bam_destroy1(b), bam_destroy1(b1);
}
