#ifndef BMF_VETTER_H
#define BMF_VETTER_H

#include <stdio.h>
#include <stdint.h>
#include <getopt.h>
#include "include/igamc_cephes.h"
#include "htslib/khash.h"
#include "htslib/vcf.h"
#include "htslib/faidx.h"
#include "htslib/sam.h"
#include "htslib/vcf.h"
#include "htslib/tbx.h"
#include "dlib/bed_util.h"
#include "dlib/bam_util.h"
#include "dlib/misc_util.h"
#include "dlib/vcf_util.h"

#include <vector>
#include <set>
#ifdef __GNUC__
#	include <parallel/algorithm>
#else
#	include <algorithm>
#endif

int vetter_main(int, char **);

KHASH_MAP_INIT_STR(names, const bam_pileup1_t *)

static int vet_func(void *data, bam1_t *b);
extern bam_plp_t bam_plp_maxcnt_init(bam_plp_auto_f func, void *data, int maxcnt);

#define DEFAULT_MAX_ALLELES 20uL

// Need to expand for new options, but I'll wait until I'm finished first.
/*
		}*/
#define VETTER_OPTIONS \
	{"min-family-agreed",		 required_argument, NULL, 'a'}, \
	{"min-family-size",		  required_argument, NULL, 's'}, \
	{"min-fraction-agreed",		 required_argument, NULL, 'f'}, \
	{"min-mapping-quality",		 required_argument, NULL, 'm'}, \
	{"min-phred-quality",		 required_argument, NULL, 'v'}, \
	{"min-count",		 required_argument, NULL, 'c'}, \
	{"min-duplex",		 required_argument, NULL, 'D'}, \
	{"min-overlap",		 required_argument, NULL, 'O'}, \
	{"out-vcf",		 required_argument, NULL, 'o'}, \
	{"bedpath",		 required_argument, NULL, 'b'}, \
	{"ref",		 required_argument, NULL, 'r'}, \
	{"padding",		 required_argument, NULL, 'p'}, \
	{"skip-secondary", no_argument, NULL, '2'},\
	{"skip-supplementary", no_argument, NULL, 'S'},\
	{"skip-qcfail", no_argument, NULL, 'q'},\
	{"skip-improper", no_argument, NULL, 'P'},\
	{"max-depth", required_argument, NULL, 'd'},\
	{"emit-bcf", no_argument, NULL, 'B'},\
	{0, 0, 0, 0}

const char *bmf_header_lines[] =  {
		"##INFO=<ID=BMF_VET,Number=A,Type=Integer,Description=\"1 if the variant passes vetting, 0 otherwise.\">",
		"##INFO=<ID=BMF_UNIOBS,Number=A,Type=Integer,Description=\"Number of unique observations supporting variant at position.\">",
		"##INFO=<ID=BMF_DUPLEX,Number=A,Type=Integer,Description=\"Number of duplex reads supporting variant at position.\">",
		"##INFO=<ID=BMF_FAIL,Number=A,Type=Integer,Description=\"Number of unique observations at position failing filters.\">",
		"##INFO=<ID=DUPLEX_DEPTH,Number=1,Type=Integer,Description=\"Number of duplex reads passing filters at position.\">",
		"##INFO=<ID=DISC_OVERLAP,Number=1,Type=Integer,Description=\"Number of read pairs at position with discordant base calls.\">",
		"##INFO=<ID=OVERLAP,Number=1,Type=Integer,Description=\"Number of overlapping read pairs combined into single observations at position.\">"
};

extern void bed_destroy(void *);
extern std::vector<khiter_t> make_sorted_keys(khash_t(bed) *h);

#endif /* BMF_VETTER_H */
