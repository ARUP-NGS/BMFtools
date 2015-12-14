#include "bam_util.h"

bam_plp_t bam_plp_maxcnt_init(bam_plp_auto_f func, void *data, int maxcnt)
{
	bam_plp_t iter = bam_plp_init(func, data);
	bam_plp_set_maxcnt(iter, maxcnt);
    return iter;
}
