#ifndef BMF_VETTER_H
#define BMF_VETTER_H

#include <stdio.h>
#include <stdint.h>
#include <stddef.h>
#include <stdio.h>
#include <omp.h>
#include <inttypes.h>
#include <getopt.h>
#include "htslib/khash.h"
#include "htslib/vcf.h"
#include "htslib/faidx.h"
#include "htslib/sam.h"
#include "htslib/vcf.h"
#include "htslib/tbx.h"
#include "bed_util.h"
#include "bam_util.h"

static int vet_func(void *data, bam1_t *b);
extern bam_plp_t bam_plp_maxcnt_init(bam_plp_auto_f func, void *data, int maxcnt);

// Setup needed for pileup engine
#ifndef BAM_PLP_DEC
#define BAM_PLP_DEC
extern void bam_plp_init_overlaps(bam_plp_t);

typedef struct {
    int k, x, y, end;
} cstate_t;

typedef struct __linkbuf_t {
    bam1_t b;
    int32_t beg, end;
    cstate_t s;
    struct __linkbuf_t *next;
} lbnode_t;

typedef struct {
    int cnt, n, max;
    lbnode_t **buf;
} mempool_t;

KHASH_MAP_INIT_STR(olap_hash, lbnode_t *)
typedef khash_t(olap_hash) olap_hash_t;

struct __bam_plp_t {
    mempool_t *mp;
    lbnode_t *head, *tail, *dummy;
    int32_t tid, pos, max_tid, max_pos;
    int is_eof, max_plp, error, maxcnt;
    uint64_t id;
    bam_pileup1_t *plp;
    // for the "auto" interface only
    bam1_t *b;
    bam_plp_auto_f func;
    void *data;
    olap_hash_t *overlaps;
};

#endif /* BAM_PLP_DEC */


#define VETTER_OPTIONS \
    {"min-family-agreed",         required_argument, NULL, 'a'}, \
    {"min-family-size",         required_argument, NULL, 's'}, \
    {"min-fraction-agreed",         required_argument, NULL, 'f'}, \
    {"min-mapping-quality",         required_argument, NULL, 'm'}, \
    {"min-phred-quality",         required_argument, NULL, 'v'}, \
    {"out-vcf",         required_argument, NULL, 'o'}, \
    {"bedpath",         required_argument, NULL, 'b'}, \
    {"ref",         required_argument, NULL, 'r'}, \
    {"padding",         required_argument, NULL, 'p'}, \
	{0, 0, 0, 0}

typedef struct vparams {
	uint32_t minFA; // Minimum Family Members Agreed on base
	uint32_t minFM; // Minimum family size
	double minFR; // Minimum  fraction agreed on base
	uint32_t minMQ; // Minimum mapping quality to include
	uint32_t minPV; // Minimum PV score to include
	uint64_t flag;
} vparams_t;

typedef struct vetplp_conf {
	bam_plp_auto_f func;
	vparams_t params;
	int n_regions;
	samFile *bam;
	bam_hdr_t *bh;
	hts_itr_t *bam_iter;
	hts_idx_t *bi; // Bam index
	hts_idx_t *vi;
	khash_t(bed) *bed; // Really khash_t(bed) *
	faidx_t *fai;
	vcfFile *vin;
	vcfFile *vout;
	bcf_hdr_t *vh;
	tbx_t *tbx;
	bam_plp_t *pileup;
	char *contig; // Holds the string for contig
	int32_t last_ref_tid; // Holds the transcript ID for the contig string.
	uint32_t minFA; // Minimum Family Members Agreed on base
	uint32_t minFM; // Minimum family size
	double minFR; // Minimum  fraction agreed on base
	uint32_t minMQ; // Minimum mapping quality to include
	uint32_t minPV; // Minimum PV score to include
	uint64_t flag;
} vetplp_conf_t;

#define SKIP_IMPROPER 4096
#define BAM_FETCH_BUFFER 150

extern void *bed_read(const char *fn);


typedef struct vetter_settings {

	char bam_path[200]; // Path to input bam
	char out_vcf_path[200]; // Path to output vcf path
	char in_vcf_path[200]; // Path to input vcf
	char bed_path[200]; // Path to bedfile
	char ref_path[200];
	char bam_rmode[4];
	char vcf_rmode[4];
	char vcf_wmode[4];
	vetplp_conf_t conf;
} vetter_settings_t;

extern void bed_destroy(void *);

#endif /* BMF_VETTER_H */
