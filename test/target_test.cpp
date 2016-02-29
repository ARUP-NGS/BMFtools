#include "htslib/sam.h"
#include "src/bmf_target.h" // Provides bed, compiler utils
#include "bam.h"

extern target_counts_t target_core(const char *bedpath, const char *bampath,
                                uint32_t padding, uint32_t minMQ, uint64_t notification_interval);
/*
#ifdef assert
#undef assert
#define assert(x) \
    do {\
        if(!(x)) {\
            LOG_EXIT((char *)"Assert failed! Assertion: '%s'.", #x);\
        }\
    } while(0)
#endif
*/

int main(int c, char **argv)
{
    target_counts_t counts = target_core((char *)"test/target_test.bed", (char *)"test/target_test.bam", 0u, 0u, 1000000);
    assert(counts.count == 2185uL);
    assert(counts.n_skipped == 0uL);
    assert(counts.target == 2032uL);
    return 0;
}
