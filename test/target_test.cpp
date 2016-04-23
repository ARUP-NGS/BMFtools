#define __STDC_LIMIT_MACROS
#include <assert.h>
#include <stdint.h>
//#include "htslib/sam.h"
#include "src/bmf_target.h"
#include "dlib/logging_util.h"

int main(int c, char **argv)
{
    BMF::target_counts_t counts = BMF::target_core((char *)"test/target_test.bed", (char *)"test/target_test.bam", 0u, 0u, 1000000);
    assert(counts.count == 2185uL);
    assert(counts.n_skipped == 0uL);
    if(counts.target != 2006uL)
        LOG_EXIT("counts.target %lu rather than expected.\n", counts.target, 2006uL);
    return 0;
}
