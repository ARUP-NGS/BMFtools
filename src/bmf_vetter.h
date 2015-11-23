#ifndef BMF_VETTER_H
#define BMF_VETTER_H

#include <stdio.h>
#include <stdint.h>
#include <stddef.h>
#include <stdio.h>
#include <omp.h>
#include <inttypes.h>
#include "htslib/vcf.h"
#include "htslib/faidx.h"

void vs_destroy(vetter_settings_t *s)
{
}

#define VETTER_OPTIONS(o1, o2, o3, o4, o5) \
    {"minFA",         required_argument, NULL, 'a'}, \
    {"minFM",         required_argument, NULL, 's'}, \
    {"minFR",         required_argument, NULL, 'p'}, \
    {"minMQ",         required_argument, NULL, 'm'}, \
    {"minPV",         required_argument, NULL, 'b'}, \
    {"ref",         required_argument, NULL, 'f'} \
	{0, 0, 0, 0}


typedef struct vetter_settings {
	int minFA; // Minimum Family Members Agreed on base
	int minFM; // Minimum family size
	int minFR; // Minimum  fraction agreed on base
	int minMQ; // Minimum mapping quality to include
	int minPV; // Minimum PV score to include
	char bam_path[200]; // Path to input vcf
	char out_vcf_path[200]; // Path to output vcf path
	char in_vcf_path[200];
	char ref_path[200];
	faidx_t *fai;
	samFile *bam;
	bam_hdr_t *bh;
	vcfFile *vin;
	vcfFile *vout;
	bcf_hdr_t *vh;
} vetter_settings_t;

#endif /* BMF_VETTER_H */
