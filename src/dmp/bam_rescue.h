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

#ifndef bam_sort_core_key
#define bam_sort_core_key(a) (uint64_t)((uint64_t)a->core.tid<<32|(a->core.pos+1)<<1|bam_is_rev(a))
#endif

#ifndef bam_sort_mate_key
#define bam_sort_mate_key(a) (uint64_t)((uint64_t)a->core.mtid<<32|a->core.mpos+1)
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
#define IS_REVERSE(bam) ((bam)->core.flag&BAM_FREVERSE != 0)
#endif

#ifndef IS_MATE_REVERSE
#define IS_MATE_REVERSE(bam) ((bam)->core.flag&BAM_FMREVERSE != 0)
#endif

#ifndef IS_READ2
#define IS_READ2(bam) ((bam)->core.flag&BAM_FREAD2 != 0)
#endif

#ifndef IS_READ1
#define IS_READ1(bam) ((bam)->core.flag&BAM_FREAD1 != 0)
#endif

#ifndef ABORT_MISSION
#define ABORT_MISSION(message) do {fprintf(stderr, message); exit(EXIT_FAILURE);} while(0)
#endif
