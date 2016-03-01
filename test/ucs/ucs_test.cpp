#include "dlib/bam_util.h"
#include <assert.h>

int main(int argc, char **argv)
{
	dlib::BamHandle handle("ucs_test.bam");
	handle.next();

    while(strcmp(bam_get_qname(handle.rec), "CACAAGTACCCATAATAA")) handle.next();

    assert(get_unclipped_start(handle.rec) == 11167320);
    return EXIT_SUCCESS;
}
