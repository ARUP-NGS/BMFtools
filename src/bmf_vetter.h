#ifndef BMF_VETTER_H
#define BMF_VETTER_H

#include <stdio.h>
#include <stdint.h>
#include <stddef.h>
#include <stdio.h>
#include <omp.h>
#include <inttypes.h>
#include "htslib/vcf.h"

void vs_destroy(vetter_settings_t *s)
{
	if(s->bh) {

	}
}


typedef struct vetter_settings {
	int minMQ; // Minimum mapping quality to include
	int minBQ; // Minimum BQ (quality string) to include
	int minPV; // Minimum PV score to include
	int minFA; // Minimum Family Members Agreed on base
	int minFM; // Minimum family size
	int minFR; // Minimum  fraction agreed on base
	char bam_path[200]; // Path to input vcf
	char out_vcf_path[200]; // Path to output vcf path
	char in_vcf_path[200];
	samFile *bam;
	bam_hdr_t *bh;
	vcfFile *vin;
	vcfFile *vout;
	bcf_hdr_t *vhin;
	bcf_hdr_t *vhout;
} vetter_settings_t;

#endif /* BMF_VETTER_H */
