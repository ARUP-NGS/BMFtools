#ifndef BMF_VETTER_H
#define BMF_VETTER_H

#include <stdio.h>
#include <stdint.h>
#include <stddef.h>
#include <stdio.h>
#include <omp.h>
#include <inttypes.h>
#include <getopt.h>
#include "htslib/vcf.h"
#include "htslib/faidx.h"
#include "htslib/sam.h"
//#include "bed_util.h"


#define VETTER_OPTIONS \
    {"min-family-agreed",         required_argument, NULL, 'a'}, \
    {"min-family-size",         required_argument, NULL, 's'}, \
    {"min-fraction-agreed",         required_argument, NULL, 'f'}, \
    {"min-mapping-quality",         required_argument, NULL, 'm'}, \
    {"min-phred-quality",         required_argument, NULL, 'p'}, \
    {"in-vcf",         required_argument, NULL, 'v'}, \
    {"out-vcf",         required_argument, NULL, 'o'}, \
    {"bedpath",         required_argument, NULL, 'b'}, \
    {"ref",         required_argument, NULL, 'r'}, \
	{0, 0, 0, 0}

typedef struct vparams {
	uint32_t minFA; // Minimum Family Members Agreed on base
	uint32_t minFM; // Minimum family size
	double minFR; // Minimum  fraction agreed on base
	uint32_t minMQ; // Minimum mapping quality to include
	uint32_t minPV; // Minimum PV score to include
} vparams_t;


typedef struct vetter_settings {

	char bam_path[200]; // Path to input bam
	char out_vcf_path[200]; // Path to output vcf path
	char in_vcf_path[200]; // Path to input vcf
	char bed_path[200]; // Path to bedfile
	char ref_path[200];
	faidx_t *fai;
	samFile *bam;
	bam_hdr_t *bh;
	vcfFile *vin;
	vcfFile *vout;
	bcf_hdr_t *vh;
	void *bed; // Really reghash_t *
	vparams_t params;
} vetter_settings_t;

extern void *bed_read(char *fn);
extern void bed_destroy(void *);

#endif /* BMF_VETTER_H */
