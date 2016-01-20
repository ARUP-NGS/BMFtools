#include "dlib/compiler_util.h"
#include "dlib/bed_util.h"
#include "htslib/sam.h"
#include "bam.h"

int main(int c, char **argv)
{
	samFile *fp = sam_open("test/target_test.bam", "r");
	bam_hdr_t *hdr = sam_hdr_read(fp);
	khash_t(bed) *bed = parse_bed_hash("test/target_test.bed", hdr, 0);
	bam1_t *b = bam_init1();
	uint64_t target = 0, count = 0, n_skipped = 0;
	while (LIKELY((c = sam_read1(fp, hdr, b)) >= 0)) {
		if((b->core.flag & (2820))) { // 2820 is unmapped, secondary, supplementary, qcfail
			++n_skipped;
			continue;
		}
		if(bed_test(b, bed)) ++target;
		if(++count % 1000000 == 0) {
			LOG_INFO("Records read: %lu.\n", count);
		}
	}
	bam_destroy1(b);
	bam_hdr_destroy(hdr);
	sam_close(fp);
	bed_destroy_hash(bed);
	//LOG_INFO("Read: %lu. Skipped: %lu. On target: %lu.\n", count, n_skipped, target);
	assert(count == 2185uL);
	assert(n_skipped == 0uL);
	assert(count == 2032L);
	return 0;
}
