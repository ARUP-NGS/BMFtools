#ifndef SEQ_UTIL_H
#define SEQ_UTIL_H
#include "compiler_util.h"
#include "khash.h"
#include "stdint.h"

KHASH_MAP_INIT_INT(shen, uint32_t)

#define nqscores 39uL
#define MAX_PV 3117 // Maximum seen with 

#define check_bam_tag(bamrec, tag) \
	do {\
		if(!bam_aux_get(bamrec, tag)) {\
		fprintf(stderr, "[E:%s] Required bam tag '%s' not found. Abort mission!\n", __func__, tag),\
		exit(EXIT_FAILURE);\
		}\
	} while(0)

#endif
