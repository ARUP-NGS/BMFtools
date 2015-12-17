#ifndef PAIR_UTIL_H
#define PAIR_UTIL_H
#include "io_util.h"
#include <stdlib.h>
#include <stdio.h>
#include <zlib.h>
#include <unistd.h>
#include "sam_opts.h"
#include "sam.h"
#include "bam.h"
#include "charcmp.h"
#include "compiler_util.h"

typedef void (*pair_fn)(bam1_t *b,bam1_t *b1);
typedef void (*single_fn)(bam1_t *b);
typedef void (*single_aux)(bam1_t *b, void *data);
typedef int (*single_aux_check)(bam1_t *b, void *data);
static inline void add_unclipped_mate_starts(bam1_t *b1, bam1_t *b2);
void abstract_pair_iter(samFile *in, bam_hdr_t *hdr, samFile *ofp, pair_fn function);
void abstract_single_filter(samFile *in, bam_hdr_t *hdr, samFile *out, single_aux_check function, void *data);
void abstract_single_data(samFile *in, bam_hdr_t *hdr, samFile *out, single_aux function, void *data);
void abstract_single_iter(samFile *in, bam_hdr_t *hdr, samFile *out, single_fn function);

enum htseq {
	HTS_A = 1,
	HTS_C = 2,
	HTS_G = 4,
	HTS_T = 8,
	HTS_N = 15
};


inline void process_mei_tag(bam1_t *b) {
	const uint8_t *tag_ptr = bam_aux_get(b, "ME");
	if(UNLIKELY(!tag_ptr)) {
		fprintf(stderr, "[E:%s] Expected ME tag not present. Abort mission! Qname: %s.", __func__,
				(char *)bam_get_qname(b));
	}
	if(bam_aux2i(tag_ptr)) {
		b->core.pos = b->core.mpos;
		b->core.tid = b->core.mtid;
		b->core.flag |= BAM_FUNMAP;
	}
	else {
		b->core.mpos = b->core.pos;
		b->core.mtid = b->core.tid;
		b->core.flag |= BAM_FMUNMAP;
	}
	return;
}

static inline void add_unclipped_mate_starts(bam1_t *b1, bam1_t *b2) {
	uint32_t i, offset1 = 0, offset2 = 0, ucs1 = b1->core.pos, ucs2 = b2->core.pos;
	const uint32_t *cigar1 = bam_get_cigar(b1);
	const uint32_t *cigar2 = bam_get_cigar(b2);
	for(i = 0; i < b1->core.n_cigar; ++i) {
		if(!(cigar1[i]&0xf)) // 'M' in cigar.
			break;
		else
			offset1 += cigar1[i] >> BAM_CIGAR_SHIFT;
	}
	for(i = 0; i < b2->core.n_cigar; ++i) {
		if(!(cigar2[i]&0xf))
			break;
		else
			offset2 += cigar2[i] >> BAM_CIGAR_SHIFT;
	}
	if(b1->core.flag & BAM_FREVERSE)
		ucs1 += offset1;
	else
		ucs1 -= offset1;
	if(b2->core.flag & BAM_FREVERSE)
		ucs2 += offset2;
	else
		ucs2 -= offset2;
	bam_aux_append(b2, "MU", 'I', sizeof(uint32_t), (uint8_t *)&ucs1);
	bam_aux_append(b1, "MU", 'I', sizeof(uint32_t), (uint8_t *)&ucs2);
	bam_aux_append(b2, "SU", 'I', sizeof(uint32_t), (uint8_t *)&ucs2);
	bam_aux_append(b1, "SU", 'I', sizeof(uint32_t), (uint8_t *)&ucs1);
	return;
}

#define check_bam_tag(bamrec, tag) \
	do {\
		if(!bam_aux_get(bamrec, tag)) {\
		fprintf(stderr, "[E:%s] Required bam tag '%s' not found. Abort mission!\n", __func__, tag),\
		exit(EXIT_FAILURE);\
		}\
	} while(0)

extern void bam_plp_set_maxcnt(bam_plp_t, int);

#endif
