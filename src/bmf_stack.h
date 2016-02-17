#ifndef BMF_STACK_H
#define BMF_STACK_H
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

KHASH_MAP_INIT_STR(names, const bam_pileup1_t *)

int stack_main(int argc, char *argv[]);


#endif /* BMF_STACK_H */
