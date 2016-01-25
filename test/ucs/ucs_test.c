#include "dlib/bam_util.h"
#include <assert.h>

int main(int argc, char **argv)
{
	samFile *fp = sam_open("ucs_test.bam", "rb");
	bam_hdr_t *hdr = sam_hdr_read(fp);
	bam1_t *b = bam_init1();

	sam_read1(fp, hdr, b);
	while(strcmp(bam_get_qname(b), "CACAAGTACCCATAATAA")) {
		sam_read1(fp, hdr, b);
	}

	int32_t ucs = get_unclipped_start(b);
	LOG_INFO("ucs: %i.\n", ucs);
	assert(ucs == 11167320);
	bam_destroy1(b);
	sam_close(fp);
	bam_hdr_destroy(hdr);
	return EXIT_SUCCESS;
}
