#ifndef BMF_CAP_H
#define BMF_CAP_H
#include "dlib/bam_util.h"
#include "include/sam_opts.h"

int cap_qscore_main(int argc, char *argv[]);

typedef struct cap_settings {
    int32_t minFM;
    uint32_t minPV;
    double minFrac;
    int cap;
    int dnd; // Do not disturb. Use uncapped unless greater than given minPV.
} cap_settings;


static inline int cap_bam_dnd(bam1_t *b, cap_settings *settings) {
    int i;
    uint32_t *PV = (uint32_t *)array_tag(b, "PV");
    uint32_t *FA = (uint32_t *)array_tag(b, "FA");
    const int FM = bam_aux2i(bam_aux_get(b, "FM"));
    if(FM < settings->minFM)
        return 1;
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
    uint32_t *const PV = (uint32_t *)array_tag(b, "PV");
    uint32_t *const FA = (uint32_t *)array_tag(b, "FA");
    const int FM = bam_aux2i(bam_aux_get(b, "FM"));
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

#endif /* BMF_CAP_H */
