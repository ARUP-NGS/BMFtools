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

namespace BMF {

	class UniqueObservation;

	class VCFPos {
		// Let's consider using VCFLib!
		std::vector<UniqueObservation> templates;
		int32_t pos;
		int32_t tid;
	};

	class UniqueObservation {
	std::string qname;
	int cycle;
	uint32_t size;
	uint32_t quality;
	uint32_t agreed;
	uint32_t mq;
	double pvalue;
	int flag; // May turn into a flag
	char base;
	public:
		UniqueObservation(const bam_pileup1_t plp):
			qname(bam_get_qname(plp.b)),
			cycle(arr_qpos(&plp)),
			size(bam_aux2i(bam_aux_get(plp.b, "FM"))),
			quality(((uint32_t *)array_tag(plp.b, "PV"))[cycle]),
			agreed(((uint32_t *)array_tag(plp.b, "FA"))[cycle]),
			mq(plp.b->core.qual),
			pvalue(std::pow(10, quality - 0.1)),
			flag(plp.b->core.flag),
			base(seq_nt16_str[bam_seqi(bam_get_seq(plp.b), plp.qpos)])
		{
		}
		void add_obs(const bam_pileup1_t plp) {
			LOG_ASSERT(strcmp(qname.c_str(), bam_get_qname(plp.b)) == 0);
			size += bam_aux2i(bam_aux_get(plp.b, "FM"));
			char _base = seq_nt16_str[bam_seqi(bam_get_seq(plp.b), plp.qpos)];
			int _cycle = arr_qpos(&plp);
			if(_base == base) {
				agreed += ((uint32_t *)array_tag(plp.b, "FM"))[_cycle];
				pvalue = igamc_pvalues(2, (double)quality + ((uint32_t *)array_tag(plp.b, "PV"))[_cycle]);
				quality = pvalue_to_phred(pvalue);
			}
		}
	};

}

KHASH_MAP_INIT_STR(names, const bam_pileup1_t *)

int stack_main(int argc, char *argv[]);


#endif /* BMF_STACK_H */
