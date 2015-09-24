#include "stdint.h"
#include "stdio.h"

typedef struct rescue_settings {
    FILE *fp;
    int write_nc2bam; // Write the reads with the number of nucleotides changed to bam rather than to a fastq file.
    int hd_thresh; // Threshold for hamming distance. If < hd_thresh, reads are considered to be from the same family.
    int is_se; // Is single end
    int force_se; // Is single end
    char *fq_fname;
} rescue_settings_t;


#ifndef REVERSE_BOOL
#define REVERSE_BOOL(bam) ((bam)->core.flag&BAM_FREVERSE != 0)
#endif

#ifndef MREVERSE_BOOL
#define MREVERSE_BOOL(bam) ((bam)->core.flag&BAM_FMREVERSE != 0)
#endif

#ifndef IS_READ2_BOOL
#define IS_READ2_BOOL(bam) ((bam)->core.flag&BAM_FREAD2 != 0)
#endif

#ifndef IS_READ1_BOOL
#define IS_READ1_BOOL(bam) ((bam)->core.flag&BAM_FREAD1 != 0)
#endif

#ifndef forever
#define forever for(;;)
#endif

#ifndef bam_is_r1
#define bam_is_r1(b) (((b)->core.flag&BAM_FREAD1) != 0)
#endif

#ifndef bam_is_r2
#define bam_is_r2(b) (((b)->core.flag&BAM_FREAD2) != 0)
#endif

#ifndef BMF_SORT_ORDER
#define BMF_SORT_ORDER 2
#endif

#ifndef SAMTOOLS_SORT_ORDER
#define SAMTOOLS_SORT_ORDER 0
#endif

#ifndef QNAME_SORT_ORDER
#define QNAME_SORT_ORDER 1
#endif

#ifndef IS_REVERSE
#define IS_REVERSE(bam) ((bam)->core.flag&BAM_FREVERSE)
#endif

#ifndef IS_MATE_REVERSE
#define IS_MATE_REVERSE(bam) ((bam)->core.flag&BAM_FMREVERSE)
#endif

#ifndef IS_READ2
#define IS_READ2(bam) ((bam)->core.flag&BAM_FREAD2)
#endif

#ifndef IS_READ1
#define IS_READ1(bam) ((bam)->core.flag&BAM_FREAD1)
#endif

typedef struct bmf_cmpkey {
    uint64_t key1;
    uint64_t key2;
} bmf_cmpkey_t;


#ifndef BMF_KEY_EQ
#define BMF_KEY_EQ(var1, var2) (var1.key1 == var2.key1 && var1.key2 == var2.key2)
#endif

#ifndef BMF_KEY_LT
#define BMF_KEY_LT(var1, var2) ((var1.key1 < var2.key1) || (var1.key2 < var2.key2))
#endif

static inline void ARR_SETKEY(bam1_t *bam, bmf_cmpkey_t key)
{
    key.key1 = (((uint64_t)bam->core.tid) >> 56 | ((uint64_t)bam->core.pos) >> 28 |
               ((uint64_t)bam->core.mtid));
    key.key2 = ((uint64_t)bam->core.mpos >> 36 | ((uint64_t)REVERSE_BOOL(bam)) >> 35 |
                ((uint64_t)MREVERSE_BOOL(bam)) >> 34 | ((uint64_t)IS_READ1_BOOL(bam)) >> 33 |
                  ((uint64_t)IS_READ1_BOOL(bam)) >> 32 | ((uint64_t)IS_READ1_BOOL(bam)) >> 31);
    return;
}

