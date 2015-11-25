#ifndef SORT_UTIL_H
#define SORT_UTIL_H

#ifndef bam_is_r1
#define bam_is_r1(b) (!!((b)->core.flag&BAM_FREAD1))
#endif

#ifndef bam_is_r2
#define bam_is_r2(b) (!!((b)->core.flag&BAM_FREAD2))
#endif


#ifndef ucs_se_sort_key
#define ucs_se_sort_key(a) (uint64_t)((uint64_t)a->core.tid<<32|(bam_aux2i(bam_aux_get(b, "SU"))+1)<<2|bam_is_rev(a))
#endif

#ifndef ucs_sort_mate_key
#define ucs_sort_mate_key(a) (uint64_t)((uint64_t)a->core.mtid<<32|(bam_aux2i(bam_aux_get(b, "MU")) + 1)<<2|(!!(a->core.flag & BAM_FMREVERSE)))
#endif

#ifndef ucs_sort_core_key
#define ucs_sort_core_key(a) (uint64_t)((uint64_t)a->core.tid<<32|(bam_aux2i(bam_aux_get(b, "SU"))+1)<<2|bam_is_rev(a)<<1|bam_is_r1(a))
#endif

#ifndef bmfsort_core_key
#define bmfsort_core_key(a) (uint64_t)((uint64_t)a->core.tid<<32|(a->core.pos+1)<<2|bam_is_rev(a)<<1|bam_is_r1(a))
#endif

#ifndef bmfsort_mate_key
#define bmfsort_mate_key(a) (uint64_t)((uint64_t)a->core.mtid<<32|(a->core.mpos+1)<<1|(!!(a->core.flag & BAM_FMREVERSE)))
#endif


#ifndef IS_REVERSE
#define IS_REVERSE(bam) ((bam)->core.flag&BAM_FREVERSE != 0)
#endif

#ifndef IS_MATE_REVERSE
#define IS_MATE_REVERSE(bam) ((bam)->core.flag&BAM_FMREVERSE != 0)
#endif

#ifndef IS_READ2
#define IS_READ2(bam) (((bam)->core.flag&BAM_FREAD2) != 0)
#endif

#ifndef IS_READ1
#define IS_READ1(bam) (((bam)->core.flag&BAM_FREAD1) != 0)
#endif

#endif /* SORT_UTIL_H */
