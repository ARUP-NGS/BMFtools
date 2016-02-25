#ifndef UNIQUE_OBS_H
#define UNIQUE_OBS_H

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
#include <unordered_map>
#include <algorithm>

namespace BMF {

	class SampleVCFPos;

	class UniqueObservation {
	friend SampleVCFPos;
	std::string qname;
	int cycle1;
	int cycle2; // Masked, from other read, if it was found.
	int discordant;
	uint32_t size;
	uint32_t quality;
	uint32_t agreed;
	uint32_t mq1;
	uint32_t mq2;
	uint32_t rv;
	int is_duplex1;
	int is_duplex2;
	int is_overlap;
	double pvalue;
	int flag; // May turn into a flag
	char base1;
	char base2; // Masked, from other read
	char base_call;
	public:
		uint32_t get_meanMQ() {
			return mq2 == (uint32_t)-1 ? mq1: ((mq2 + mq1 + 0.5) / 2);
		}
		double get_FA() {
			return (double)agreed / size;
		}
		int get_overlap() {return is_overlap;}
		int get_duplex() {
			return is_duplex1 + (is_duplex2 != -1 ? is_duplex2: 0);
		}
		UniqueObservation() {
			memset(this, 0, sizeof(*this));
		}
		inline UniqueObservation(const bam_pileup1_t& plp):
			qname(bam_get_qname(plp.b)),
			cycle1(arr_qpos(&plp)),
			cycle2(-1),
			discordant(-1),
			size(bam_aux2i(bam_aux_get(plp.b, "FM"))),
			quality(((uint32_t *)array_tag(plp.b, "PV"))[cycle1]),
			agreed(((uint32_t *)array_tag(plp.b, "FA"))[cycle1]),
			mq1(plp.b->core.qual),
			mq2((uint32_t)-1),
			rv((uint32_t)bam_aux2i(bam_aux_get(plp.b, "RV"))),
			is_duplex1(bam_aux2i(bam_aux_get(plp.b, "DR"))),
			is_duplex2(-1),
			is_overlap(0),
			pvalue(std::pow(10, quality - 0.1)),
			flag(plp.b->core.flag),
			base1(seq_nt16_str[bam_seqi(bam_get_seq(plp.b), plp.qpos)]),
			base2('\0'),
			base_call(base1)
		{
		}
		void add_obs(const bam_pileup1_t& plp);
	};

	class PairVCFPos;

	class SampleVCFPos {
		friend PairVCFPos;
		std::unordered_map<char, std::vector<UniqueObservation *>> templates;
		size_t size;
		int32_t pos;
		int32_t tid;
	public:
		SampleVCFPos(std::unordered_map<char *, UniqueObservation> obs, int32_t _tid, int32_t _pos);
		void to_bcf(bcf1_t *vrec, bcf_hdr_t *hdr, char refbase);
	};

	class PairVCFPos {
		SampleVCFPos tumor;
		SampleVCFPos normal;
	public:
		void to_bcf(bcf1_t *vrec, bcf_hdr_t *hdr, char refbase);
		PairVCFPos(std::unordered_map<char *, UniqueObservation> tobs, std::unordered_map<char *, UniqueObservation> nobs,
					int32_t _tid, int32_t _pos):
						tumor(tobs, _tid, _pos),
						normal(nobs, _tid, _pos) {
		}
	};
}

#endif
