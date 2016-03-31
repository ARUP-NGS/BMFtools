#ifndef BMF_VETTER_H
#define BMF_VETTER_H

#include "dlib/bam_util.h"
#include "dlib/bed_util.h"
#include "dlib/misc_util.h"
#include "dlib/vcf_util.h"
#include <stdio.h>
#include <stdint.h>
#include <getopt.h>
#include "include/igamc_cephes.h"
#include "include/bam.h"
#include "htslib/khash.h"
#include "htslib/vcf.h"
#include "htslib/faidx.h"
#include "htslib/sam.h"
#include "htslib/vcf.h"
#include "htslib/tbx.h"

#include <vector>
#include <set>
#include <algorithm>

namespace BMF {

    int vet_main(int, char **);

    KHASH_MAP_INIT_STR(names, const bam_pileup1_t *)

    uint64_t NUM_PREALLOCATED_ALLELES = 4uL;
    const char *bmf_header_lines[] =  {
            "##INFO=<ID=BMF_VET,Number=A,Type=Integer,Description=\"1 if the variant passes vetting, 0 otherwise.\">",
            "##INFO=<ID=BMF_UNIOBS,Number=A,Type=Integer,Description=\"Number of unique observations supporting variant at position.\">",
            "##INFO=<ID=BMF_DUPLEX,Number=A,Type=Integer,Description=\"Number of duplex reads supporting variant at position.\">",
            "##INFO=<ID=BMF_FAIL,Number=A,Type=Integer,Description=\"Number of unique observations at position failing filters.\">",
            "##INFO=<ID=DUPLEX_DEPTH,Number=1,Type=Integer,Description=\"Number of duplex reads passing filters at position.\">",
            "##INFO=<ID=DISC_OVERLAP,Number=1,Type=Integer,Description=\"Number of read pairs at position with discordant base calls.\">",
            "##INFO=<ID=OVERLAP,Number=1,Type=Integer,Description=\"Number of overlapping read pairs combined into single observations at position.\">"
    };

}

//extern std::vector<khiter_t> dlib::make_sorted_keys(khash_t(bed) *h);

#endif /* BMF_VETTER_H */
