#ifndef CAP_QSCORE_H
#define CAP_QSCORE_H
#include "bam_util.h"

int cap_qscore_main(int argc, char *argv[]);

typedef struct cap_settings {
	uint32_t minFM;
	uint32_t minPV;
	double minFrac;
	char cap;
	int dnd; // Do not disturb. Use uncapped unless greater than given minPV.
} cap_settings;


static inline int cap_bam_dnd(bam1_t *b, cap_settings *settings) {
    int i;
	uint32_t *PV = (uint32_t *)array_tag(b, (char *)"PV");
	uint32_t *FA = (uint32_t *)array_tag(b, (char *)"FA");
	const int FM = bam_aux2i(bam_aux_get(b, (char *)"FM"));
	if((uint32_t)FM < settings->minFM) {
		return 1;
	}
	const int l_qseq = b->core.l_qseq;
	char *qual = (char *)bam_get_qual(b);
	if(b->core.flag & BAM_FREVERSE)
		for(i = 0; i < l_qseq; ++i)
			if(PV[i] >= settings->minPV && (double)FA[i] / FM >= settings->minFrac)
				qual[l_qseq - 1 - i] = settings->cap;
	else
		for(i = 0; i < l_qseq; ++i)
			if(PV[i] >= settings->minPV && (double)FA[i] / FM >= settings->minFrac)
				qual[i] = settings->cap;
	return 0;
}

static inline int cap_bam_q(bam1_t *b, cap_settings *settings) {
    uint32_t *const PV = (uint32_t *)array_tag(b, (char *)"PV");
	uint32_t *const FA = (uint32_t *)array_tag(b, (char *)"FA");
	const int FM = bam_aux2i(bam_aux_get(b, (char *)"FM"));
	if(FM < (int)settings->minFM)
		return 1;
	char *qual = (char *)bam_get_qual(b);
	int i;
	const int l_qseq = b->core.l_qseq;
	if(b->core.flag & BAM_FREVERSE)
		for(i = 0; i < l_qseq; ++i)
			qual[l_qseq - 1 - i] = (PV[i] >= settings->minPV && (double)FA[i] / FM >= settings->minFrac) ? settings->cap: 2;
	else
		for(i = 0; i < l_qseq; ++i)
			qual[i] = (PV[i] >= settings->minPV && (double)FA[i] / FM >= settings->minFrac) ? settings->cap: 2;
	return 0;
}

#endif /* CAP_QSCORE_H */
