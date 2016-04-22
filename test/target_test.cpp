#define __STDC_LIMIT_MACROS
#include <assert.h>
#include <stdint.h>
//#include "htslib/sam.h"
#include "src/bmf_target.h"

int main(int c, char **argv)
{
    BMF::target_counts_t counts = BMF::target_core((char *)"test/target_test.bed", (char *)"test/target_test.bam", 0u, 0u, 1000000);
    assert(counts.count == 2185uL);
    assert(counts.n_skipped == 0uL);
    assert(counts.target == 2032uL);
    return 0;
}
