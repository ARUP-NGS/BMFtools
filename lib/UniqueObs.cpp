#include "UniqueObs.h"

namespace BMF {
	/*
	 * TODO:
	 * 1. Decide which info/format fields for each variant.
	 * 2. Decide to have a number of fields
	 */
	VCFPos::VCFPos(std::unordered_map<std::string, UniqueObservation> obs, int32_t _tid, int32_t _pos, size_t _size):
	pos(_pos),
	tid(_tid),
	size(_size){
		templates.reserve(5);
		for(auto& pair: obs) templates[pair.second.base_call].push_back(&pair.second);
	}
	void VCFPos::to_bcf(bcf1_t *vrec, bcf_hdr_t *hdr, char refbase) {
		size_t n_alleles = templates.size();
		vrec->rid = tid;
		vrec->pos = pos;
		vrec->qual = 0;
		kstring_t allele_str = {0, 0, NULL};
		ks_resize(&allele_str, 20uL);
		kputc(refbase, &allele_str);
		auto match = templates.find(refbase);
		if(match == templates.end()) {
			// No reference calls? Add appropriate info fields.
		}
		for(auto& pair: templates) {
			if(pair.first != refbase)
				 kputc(',', &allele_str), kputc((int)pair.first, &allele_str);
			// Update records
		}
		bcf_update_alleles_str(hdr, vrec, allele_str.s);
	}

	void UniqueObservation::add_obs(const bam_pileup1_t& plp) {
		LOG_ASSERT(strcmp(qname.c_str(), bam_get_qname(plp.b)) == 0);
		size += bam_aux2i(bam_aux_get(plp.b, "FM"));
		base2 = seq_nt16_str[bam_seqi(bam_get_seq(plp.b), plp.qpos)];
		cycle2 = arr_qpos(&plp);
		mq2 = (uint32_t)plp.b->core.qual;
		if(base2 == base1) {
			is_overlap = 1;
			discordant = 0;
			agreed += ((uint32_t *)array_tag(plp.b, "FM"))[cycle2];
			pvalue = agreed_pvalues(quality, ((uint32_t *)array_tag(plp.b, "PV"))[cycle2]);
			quality = pvalue_to_phred(pvalue);
			rv += bam_aux2i(bam_aux_get(plp.b, "FM"));
		} else if(base1 == 'N') {
			discordant = 0;
			// Basically, throw out
			base_call = base2;
			agreed = ((uint32_t *)array_tag(plp.b, "FM"))[cycle2];
			quality = ((uint32_t *)array_tag(plp.b, "PV"))[cycle2];
			pvalue = -10 * std::pow(10, -0.1 * quality);
			rv = (uint32_t)bam_aux2i(bam_aux_get(plp.b, "RV"));
		} else if(base2 != 'N') {
			discordant = 1;
			base_call = 'N';
			agreed = 0;
			quality = 0;
			pvalue = 1.;
		}
	}
}
