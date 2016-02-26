#include "stack.h"

namespace BMF {
	/*
	 * TODO:
	 * 1. Decide which info/format fields for each variant.
	 * 2. Decide to have a number of fields
	 */
	SampleVCFPos::SampleVCFPos(std::unordered_map<char *, UniqueObservation> obs, int32_t _tid, int32_t _pos):
	size(obs.size()),
	pos(_pos),
	tid(_tid) {
		templates.reserve(5);
		for(auto& pair: obs) templates[pair.second.base_call].push_back(&pair.second);
	}
	void SampleVCFPos::to_bcf(bcf1_t *vrec, bcf_hdr_t *hdr, char refbase) {
		std::vector<int> counts;
		std::vector<int> duplex_counts;
		std::vector<int> overlap_counts;
		int count_index = 0;
		vrec->rid = tid;
		vrec->pos = pos;
		vrec->qual = 0;
		kstring_t allele_str = {0, 0, NULL};
		ks_resize(&allele_str, 20uL);
		kputc(refbase, &allele_str);
		auto match = templates.find(refbase);
		if(match == templates.end()) {
			counts = std::vector<int>(templates.size() + 1, 0);
			duplex_counts = std::vector<int>(templates.size() + 1, 0);
			overlap_counts = std::vector<int>(templates.size() + 1, 0);
			bcf_update_info_flag(hdr, vrec, "NOREF", NULL, 1);
			// No reference calls? Add appropriate info fields.
		} else {
			counts = std::vector<int>(templates.size());
			duplex_counts = std::vector<int>(templates.size());
			overlap_counts = std::vector<int>(templates.size());
			counts[0] = match->second.size();
			for(auto uni: templates[refbase]) {
				duplex_counts[0] += uni->get_duplex();
				overlap_counts[0] += uni->get_overlap();
			}
		}
		for(auto& pair: templates) {
			if(pair.first != refbase)
				 kputc(',', &allele_str), kputc((int)pair.first, &allele_str);
			counts[++count_index] = pair.second.size();
			for(auto uni: templates[pair.first]) {
				duplex_counts[count_index] += uni->get_duplex();
				overlap_counts[count_index] += uni->get_overlap();
			}
			// Update records
		}
		bcf_update_alleles_str(hdr, vrec, allele_str.s);
		bcf_update_format(hdr, vrec, "ADP", (const void *)&counts[0], count_index, 'i');
		bcf_update_format(hdr, vrec, "ADPL", (const void *)&duplex_counts[0], count_index, 'i');
		bcf_update_format(hdr, vrec, "ADPO", (const void *)&overlap_counts[0], count_index, 'i');
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
	void PairVCFPos::to_bcf(bcf1_t *vrec, stack_aux_t *aux, int ttid, int tpos) {
		const char refbase = aux->get_ref_base(ttid, tpos);
		std::vector<char> base_calls;
		std::unordered_set<char> base_set = std::unordered_set<char>({refbase});
		for(auto& pair: tumor.templates)
			base_set.insert(pair.first);
		for(auto& pair: normal.templates)
			base_set.insert(pair.first);
		base_calls = std::vector<char>(base_set.begin(), base_set.end());
		std::sort(base_calls.begin(), base_calls.end(), [refbase](const char a, const char b) {
			return (a == refbase) ? true : (b == refbase) ? false: a < b;
		});
		std::vector<int> counts = std::vector<int>(base_calls.size() * 2);
		std::vector<int> duplex_counts = std::vector<int>(base_calls.size() * 2);
		std::vector<int> overlap_counts = std::vector<int>(base_calls.size() * 2);
		std::vector<int> allele_passes = std::vector<int>(base_calls.size() * 2);
		vrec->rid = tumor.tid;
		vrec->pos = tumor.pos;
		vrec->qual = 0;
		auto match = tumor.templates.find(refbase);
		if(match != tumor.templates.end()) {
			// Already 0-initialized if not found.
			counts[0] = match->second.size();
			for(auto uni: tumor.templates[refbase]) {
				duplex_counts[0] += uni->get_duplex();
				overlap_counts[0] += uni->get_overlap();
			}
		}
		if((match = normal.templates.find(refbase)) != normal.templates.end()) {
			counts[base_calls.size()] = match->second.size();
			for(auto uni: normal.templates[refbase]) {
				duplex_counts[base_calls.size()] += uni->get_duplex();
				overlap_counts[base_calls.size()] += uni->get_overlap();
			}
		}
		kstring_t allele_str = {0, 0, NULL};
		ks_resize(&allele_str, 8uL);
		kputc(refbase, &allele_str);
		for(unsigned i = 1; i < base_calls.size(); ++i) {
			kputc(',', &allele_str), kputc(base_calls[i], &allele_str);
			if((match = tumor.templates.find(base_calls[i])) != tumor.templates.end()) {
				counts[i] = match->second.size();
				for(auto uni: tumor.templates[base_calls[i]]) {
					duplex_counts[i] += uni->get_duplex();
					overlap_counts[i] += uni->get_overlap();
				}
			}
			if((match = normal.templates.find(base_calls[i])) != normal.templates.end()) {
				counts[i + base_calls.size()] = match->second.size();
				for(auto uni: normal.templates[base_calls[i]]) {
					duplex_counts[i + base_calls.size()] += uni->get_duplex();
					overlap_counts[i + base_calls.size()] += uni->get_overlap();
				}
			}
		}
		bcf_update_alleles_str(aux->vcf->vh, vrec, allele_str.s);
		bcf_update_format(aux->vcf->vh, vrec, "ADP", (const void *)&counts[0], counts.size(), 'i');
		bcf_update_format(aux->vcf->vh, vrec, "ADPL", (const void *)&duplex_counts[0], counts.size(), 'i');
		bcf_update_format(aux->vcf->vh, vrec, "ADPO", (const void *)&overlap_counts[0], counts.size(), 'i');
		free(allele_str.s);
	}
}