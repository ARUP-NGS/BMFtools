#ifndef BMF_VETTER_H
#define BMF_VETTER_H

#include <stdio.h>
#include <stdint.h>
#include <getopt.h>
#include "htslib/khash.h"
#include "htslib/vcf.h"
#include "htslib/faidx.h"
#include "htslib/sam.h"
#include "htslib/vcf.h"
#include "htslib/tbx.h"
#include "dlib/bed_util.h"
#include "dlib/bam_util.h"
#include "dlib/misc_util.h"

KHASH_MAP_INIT_STR(names, bam1_t *)

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

// Need to expand for new optios, but I'll wait until I'm finished first.
#define VETTER_OPTIONS \
	{"min-family-agreed",		 required_argument, NULL, 'a'}, \
	{"min-family-size",		 required_argument, NULL, 's'}, \
	{"min-fraction-agreed",		 required_argument, NULL, 'f'}, \
	{"min-mapping-quality",		 required_argument, NULL, 'm'}, \
	{"min-phred-quality",		 required_argument, NULL, 'v'}, \
	{"out-vcf",		 required_argument, NULL, 'o'}, \
	{"bedpath",		 required_argument, NULL, 'b'}, \
	{"ref",		 required_argument, NULL, 'r'}, \
	{"padding",		 required_argument, NULL, 'p'}, \
	{0, 0, 0, 0}

const char *bmf_header_lines[] =  {
		"##FORMAT=<ID=BMF_VET,Number=A,Type=Integer,Description=\"1 if the variant passes vetting, 0 otherwise.\">"
};

typedef struct {
	samFile *fp;
	hts_itr_t *iter;
	bam_hdr_t *header;
	vcfFile *vcf_fp;
	vcfFile *vcf_ofp;
	bcf_hdr_t *vcf_header;
	khash_t(bed) *bed;
	float minFR; // Minimum fraction of family members agreed on base
	float minAF; // Minimum aligned fraction
	int max_depth;
	int minFM;
	int minFA;
	int minPV;
	int minMQ;
	int minCount;
	int minDuplex;
	int minOverlap;
	int skip_improper;
	uint32_t skip_flag; // Skip reads with any bits set to true
} aux_t;

typedef struct vetplp_conf {
	samFile *bam;
	bam_hdr_t *bh;
	khash_t(bed) *bed;
	vcfFile *vin;
	vcfFile *vout;
	bcf_hdr_t *vh;
} vetplp_conf_t;

#define SKIP_IMPROPER 4096
#define BAM_FETCH_BUFFER 150

extern void *bed_read(const char *fn);


typedef struct vetter_settings {

	char *bam_path; // Path to input bam
	char *outvcf; // Path to output vcf path
	char *invcf; // Path to input vcf
	char *bed; // Path to bedfile
	char *ref_path;
	char vcf_wmode[4];
	vetplp_conf_t conf;
} vetter_settings_t;

extern void bed_destroy(void *);

#endif /* BMF_VETTER_H */
