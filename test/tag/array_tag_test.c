#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "htslib/sam.h"
#include "dlib/bam_util.h"

int main() {
	static const uint32_t arr[] = {2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,32,38,38,32,38,38,38,40,32,27,38,38,38,38,38,32,38,38,38,38,40,38,38,38,38,32,38,38,38,38,38,40,13,38,38,38,32,38,32,38,38,40,40,32,38,38,38,38,32,38,38,38,38,38,40,40,32,38,38,40,40,40,38,32,32,13,32,38,38,38,32,38,40,40,40,38,38,38,38,32,38,38,38,38,38,32,32,38,32,38,38,32,38,38,38,38,40,38,13,27,38,27,32,32,38,38,38,38,32,38,13,38,27,32,38,32,32,38,38,38,2};
	samFile *fp = sam_open("test/tag/tag_test.bam", "r");
	samFile *ofp = sam_open("-", "w");
	bam_hdr_t *h = sam_hdr_read(fp);
	if(!fp || !h)  {
		fprintf(stderr, "Could not open tag_test.bam. Abort!\n");
		exit(EXIT_FAILURE);
	}
	bam1_t *b = bam_init1();
	sam_read1(fp, h, b);
	uint32_t *PV = array_tag(b, "PV");
	assert(PV);
	assert(PV[-1] == 137);
	sam_write1(ofp, h, b);
	for(int i = 0; i < b->core.l_qseq; ++i) {
		fprintf(stderr, ",%u", PV[i]);
		assert(PV[i] = arr[i]);
	}
	bam_destroy1(b);
	bam_hdr_destroy(h);
	sam_close(fp);
	sam_close(ofp);
	return EXIT_SUCCESS;
}
