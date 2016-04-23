#include "stack.h"

#include <stdint.h>
#include <algorithm>
#include "include/igamc_cephes.h"
#include "dlib/misc_util.h"

namespace BMF {
    /*
     * TODO:
     * 1. Decide which info/format fields for each variant.
     * 2. Decide to have a number of fields
     */
    SampleVCFPos::SampleVCFPos(std::unordered_map<std::string, UniqueObservation>& obs, int32_t tid, int32_t pos):
    size(obs.size()),
    pos(pos),
    tid(tid) {
        templates.reserve(4);
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
        kstring_t allele_str = {0, 10uL, (char *)malloc(10uL * sizeof(char))};
        kputc(refbase, &allele_str);
        auto match = templates.find(refbase);
        if(match == templates.end()) {
            counts.resize(templates.size() + 1);
            duplex_counts.resize(templates.size() + 1);
            overlap_counts.resize(templates.size() + 1);
            bcf_update_info_flag(hdr, vrec, "NOREF", nullptr, 1);
            // No reference calls? Add appropriate info fields.
        } else {
            counts.resize(templates.size());
            duplex_counts.resize(templates.size());
            overlap_counts.resize(templates.size());
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
        bcf_update_format(hdr, vrec, "ADP", (const void *)counts.data(), counts.size(), 'i');
        bcf_update_format(hdr, vrec, "ADPD", (const void *)duplex_counts.data(), duplex_counts.size(), 'i');
        bcf_update_format(hdr, vrec, "ADPO", (const void *)overlap_counts.data(), overlap_counts.size(), 'i');
    }

    void UniqueObservation::add_obs(const bam_pileup1_t& plp) {
        LOG_ASSERT(strcmp(qname.c_str(), bam_get_qname(plp.b)) == 0);
        size += bam_itag(plp.b, "FM");
        base2 = seq_nt16_str[bam_seqi(bam_get_seq(plp.b), plp.qpos)];
        cycle2 = dlib::arr_qpos(&plp);
        mq2 = (uint32_t)plp.b->core.qual;
        is_reverse2 = bam_is_rev(plp.b);
        is_overlap = 1;
        if(base2 == base1) {
            discordant = 0;
            agreed += ((uint32_t *)dlib::array_tag(plp.b, "FA"))[cycle2];
            quality = agreed_pvalues(quality, ((uint32_t *)dlib::array_tag(plp.b, "PV"))[cycle2]);
            pvalue = std::pow(10, -0.1 * quality);
            rv += bam_itag(plp.b, "RV");
        } else if(base1 == 'N') {
            discordant = 0;
            // Basically, throw out
            base_call = base2;
            agreed = ((uint32_t *)dlib::array_tag(plp.b, "FA"))[cycle2];
            quality = ((uint32_t *)dlib::array_tag(plp.b, "PV"))[cycle2];
            pvalue = std::pow(10, -0.1 * quality);
            rv = (uint32_t)bam_itag(plp.b, "RV");
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
        int ambig[2] = {0, 0};
        std::unordered_set<char> base_set{refbase};
        for(auto& pair: tumor.templates) {
            if(pair.first == 'N') {
               ambig[0] = pair.second.size();
               continue;
            }
            base_set.insert(pair.first);
        }
        for(auto& pair: normal.templates) {
            if(pair.first == 'N') {
               ambig[1] = pair.second.size();
               continue;
            }
            base_set.insert(pair.first);
        }
        std::vector<char> base_calls(base_set.begin(), base_set.end());
        const size_t n_base_calls = base_calls.size();
        // Sort lexicographically AFTER putting the reference base first.
        std::sort(base_calls.begin(), base_calls.end(), [refbase](const char a, const char b) {
            return (a == refbase) ? true : (b == refbase) ? false: a < b;
        });
        std::vector<int> counts(n_base_calls * 2);
        std::vector<int> duplex_counts(n_base_calls * 2);
        std::vector<int> overlap_counts(n_base_calls * 2);
        std::vector<int> reverse_counts(n_base_calls * 2);
        std::vector<int> failed_counts(n_base_calls * 2);
        std::vector<int> allele_passes(n_base_calls * 2);
        std::vector<int> qscore_sums(n_base_calls * 2);
        std::vector<int> somatic(n_base_calls);
        vrec->rid = tumor.tid;
        vrec->pos = tumor.pos;
        vrec->qual = 0;
        vrec->n_sample = 2;
        auto match = tumor.templates.find(refbase);
        if(match != tumor.templates.end()) {
            // Already 0-initialized if not found.
            counts[0] = match->second.size();
            // auto without "reference" because auto is pointers.
            for(auto uni: tumor.templates[refbase]) {
                if(uni->get_quality() < aux->conf.minPV || uni->get_agreed() < aux->conf.minFA
                        || (float)uni->get_agreed() / uni->get_size() < aux->conf.minFR) {
                    uni->set_pass(0);
                    ++failed_counts[0];
                    continue;
                }
                duplex_counts[0] += uni->get_duplex();
                overlap_counts[0] += uni->get_overlap();
                reverse_counts[0] += uni->get_reverse();
                qscore_sums[0] += uni->get_quality();
                allele_passes[0] = (duplex_counts[0] >= aux->conf.minDuplex &&
                                    counts[0] >= aux->conf.minCount &&
                                    overlap_counts[0] >= aux->conf.minOverlap);
            }
        }
        if((match = normal.templates.find(refbase)) != normal.templates.end()) {
            counts[n_base_calls] = match->second.size();
            // auto without "reference" because auto is pointers.
            for(auto uni: normal.templates[refbase]) {
                if(uni->get_quality() < aux->conf.minPV || uni->get_agreed() < aux->conf.minFA
                       || (float)uni->get_agreed() / uni->get_size() < aux->conf.minFR) {
                    uni->set_pass(0);
                    ++failed_counts[n_base_calls];
                    continue;
                }
                duplex_counts[n_base_calls] += uni->get_duplex();
                overlap_counts[n_base_calls] += uni->get_overlap();
                reverse_counts[n_base_calls] += uni->get_reverse();
                qscore_sums[n_base_calls] += uni->get_quality();
            }
            allele_passes[n_base_calls] = (duplex_counts[n_base_calls] >= aux->conf.minDuplex &&
                                                counts[n_base_calls] >= aux->conf.minCount &&
                                                overlap_counts[n_base_calls] >= aux->conf.minOverlap);
        }
        kstring_t allele_str = {0, 0, nullptr};
        ks_resize(&allele_str, 8uL);
        kputc(refbase, &allele_str);
        for(unsigned i = 1; i < n_base_calls; ++i) {
            kputc(',', &allele_str), kputc(base_calls[i], &allele_str);
            if((match = normal.templates.find(base_calls[i])) != normal.templates.end()) {
                counts[i + n_base_calls] = match->second.size();
                for(auto uni: normal.templates[base_calls[i]]) {
                    if(uni->get_quality() < aux->conf.minPV || uni->get_agreed() < aux->conf.minFA
                            || (float)uni->get_agreed() / uni->get_size() < aux->conf.minFR) {
                        uni->set_pass(0);
                        ++failed_counts[i];
                        continue;
                    }
                    duplex_counts[i + n_base_calls] += uni->get_duplex();
                    overlap_counts[i + n_base_calls] += uni->get_overlap();
                    reverse_counts[i + n_base_calls] += uni->get_reverse();
                    qscore_sums[i + n_base_calls] += uni->get_quality();
                }
                allele_passes[i + n_base_calls] = (duplex_counts[i + n_base_calls] >= aux->conf.minDuplex &&
                                                        counts[i + n_base_calls] >= aux->conf.minCount &&
                                                        overlap_counts[i + n_base_calls] >= aux->conf.minOverlap);
            }
            if((match = tumor.templates.find(base_calls[i])) != tumor.templates.end()) {
                counts[i] = match->second.size();
                for(auto uni: tumor.templates[base_calls[i]]) {
                    if(uni->get_quality() < aux->conf.minPV || uni->get_agreed() < aux->conf.minFA
                            || (float)uni->get_agreed() / uni->get_size() < aux->conf.minFR) {
                        uni->set_pass(0);
                        ++failed_counts[i];
                        continue;
                    }
                    duplex_counts[i] += uni->get_duplex();
                    overlap_counts[i] += uni->get_overlap();
                    reverse_counts[i] += uni->get_reverse();
                    qscore_sums[i] += uni->get_quality();
                }
                allele_passes[i] = (duplex_counts[i] >= aux->conf.minDuplex &&
                                                        counts[i] >= aux->conf.minCount &&
                                                        overlap_counts[i] >= aux->conf.minOverlap);
                if(allele_passes[i] && !allele_passes[i + n_base_calls]) {
                    bcf_update_info_flag(aux->vcf.vh, vrec, "SOMATIC", nullptr, 1);
                    somatic[i] = 1;
                }
            }
        }
        const int total_depth_tumor(std::accumulate(counts.begin(), counts.begin() + n_base_calls, 0));
        const int total_depth_normal(std::accumulate(counts.begin() + n_base_calls, counts.end(), 0));
        //LOG_DEBUG("Got total depths %i,%i.\n", total_depth_tumor, total_depth_normal);
        std::vector<float> rv_fractions(reverse_counts.size());
        std::vector<float> allele_fractions(reverse_counts.size());
        for(unsigned i = 0; i < n_base_calls; ++i) {
            rv_fractions[i] = (float)counts[i] / reverse_counts[i];
            rv_fractions[i + n_base_calls] = (float)counts[i + n_base_calls] / reverse_counts[i + n_base_calls];
            allele_fractions[i] = (float)counts[i] / total_depth_tumor;
            allele_fractions[i + n_base_calls] = (float)counts[i + n_base_calls] / total_depth_normal;
        }
        bcf_update_alleles_str(aux->vcf.vh, vrec, allele_str.s), free(allele_str.s);
        bcf_update_format_int32(aux->vcf.vh, vrec, "ADP", static_cast<const void *>(counts.data()), counts.size());
        bcf_update_format_int32(aux->vcf.vh, vrec, "ADPD", static_cast<const void *>(duplex_counts.data()), duplex_counts.size());
        bcf_update_format_int32(aux->vcf.vh, vrec, "ADPO", static_cast<const void *>(overlap_counts.data()), overlap_counts.size());
        bcf_update_format_int32(aux->vcf.vh, vrec, "ADPR", static_cast<const void *>(reverse_counts.data()), reverse_counts.size());
#if !NDEBUG
        for(auto i: reverse_counts) LOG_DEBUG("Reverse count: %i.\n", i);
#endif
        bcf_update_format_float(aux->vcf.vh, vrec, "RVF", static_cast<const void *>(rv_fractions.data()), rv_fractions.size());
        bcf_update_format_int32(aux->vcf.vh, vrec, "BMF_PASS", static_cast<const void *>(allele_passes.data()), allele_passes.size());
        bcf_update_format_int32(aux->vcf.vh, vrec, "QSS", static_cast<const void *>(qscore_sums.data()), qscore_sums.size());
        bcf_update_format_int32(aux->vcf.vh, vrec, "AMBIG", static_cast<const void *>(ambig), sizeof(ambig));
        bcf_update_format_float(aux->vcf.vh, vrec, "AFR", static_cast<const void *>(allele_fractions.data()), allele_fractions.size());
        bcf_update_info_int32(aux->vcf.vh, vrec, "SOMATIC_CALL", static_cast<const void *>(somatic.data()), somatic.size());
    } /* PairVCFLine::to_bcf */

void add_hdr_lines(bcf_hdr_t *hdr, const char *lines[], size_t n) {
    while(n)
        if(bcf_hdr_append(hdr, lines[--n]))
            LOG_EXIT("Could not add header line %s. Abort!\n");
}

void add_stack_lines(bcf_hdr_t *hdr) {
    add_hdr_lines(hdr, stack_vcf_lines, sizeof(stack_vcf_lines) / sizeof(stack_vcf_lines[0]));
}

} /* namespace BMF */
