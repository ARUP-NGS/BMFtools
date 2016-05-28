#include "stack.h"

#include <stdint.h>
#include <algorithm>
#include "include/igamc_cephes.h"
#include "dlib/misc_util.h"

namespace bmf {
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
        // TODO: Re-do the paired-end mode for single-end.
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
        // TODO: Add BMF_QUANT tag to PairVCFPos.
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
        std::vector<std::vector<uint32_t>> tconfident_phreds;
        std::vector<std::vector<uint32_t>> tsuspect_phreds;
        std::vector<std::vector<uint32_t>> nconfident_phreds;
        std::vector<std::vector<uint32_t>> nsuspect_phreds;
        tconfident_phreds.reserve(n_base_calls);
        tsuspect_phreds.reserve(n_base_calls);
        nconfident_phreds.reserve(n_base_calls);
        nsuspect_phreds.reserve(n_base_calls);
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
        tsuspect_phreds.emplace_back();
        tconfident_phreds.emplace_back();
        nsuspect_phreds.emplace_back();
        nconfident_phreds.emplace_back();
        if(match != tumor.templates.end()) {
            // Already 0-initialized if not found.
            counts[0] = match->second.size();
            // auto without "reference" because auto is pointers.
            for(auto uni: tumor.templates[refbase]) {
                if(uni->get_quality() < aux->conf.minPV || uni->get_agreed() < aux->conf.minFA
                        || (float)uni->get_agreed() / uni->get_size() < aux->conf.min_fr) {
                    tsuspect_phreds[0].push_back(uni->get_quality());
                    uni->set_pass(0);
                    ++failed_counts[0];
                }
                else {
                    tconfident_phreds[0].push_back(uni->get_quality());
                    duplex_counts[0] += uni->get_duplex();
                    overlap_counts[0] += uni->get_overlap();
                    reverse_counts[0] += uni->get_reverse();
                    qscore_sums[0] += uni->get_quality();
                    allele_passes[0] = (duplex_counts[0] >= aux->conf.min_duplex &&
                                        counts[0] >= aux->conf.min_count &&
                                        overlap_counts[0] >= aux->conf.min_overlap);
                }
            }
        }
        if((match = normal.templates.find(refbase)) != normal.templates.end()) {
            counts[n_base_calls] = match->second.size();
            // auto without "reference" because auto is pointers.
            for(auto uni: normal.templates[refbase]) {
                if(uni->get_quality() < aux->conf.minPV || uni->get_agreed() < aux->conf.minFA
                       || (float)uni->get_agreed() / uni->get_size() < aux->conf.min_fr) {
                    nsuspect_phreds[0].push_back(uni->get_quality());
                    uni->set_pass(0);
                    ++failed_counts[n_base_calls];
                    continue;
                } else {
                    nconfident_phreds[0].push_back(uni->get_quality());
                    duplex_counts[n_base_calls] += uni->get_duplex();
                    overlap_counts[n_base_calls] += uni->get_overlap();
                    reverse_counts[n_base_calls] += uni->get_reverse();
                    qscore_sums[n_base_calls] += uni->get_quality();
                }
            }
            allele_passes[n_base_calls] = (duplex_counts[n_base_calls] >= aux->conf.min_duplex &&
                                           counts[n_base_calls] >= aux->conf.min_count &&
                                           overlap_counts[n_base_calls] >= aux->conf.min_overlap);
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
                            || (float)uni->get_agreed() / uni->get_size() < aux->conf.min_fr) {
                        nsuspect_phreds[i].push_back(uni->get_quality());
                        uni->set_pass(0);
                        ++failed_counts[i];
                    } else {
                        nconfident_phreds[i].push_back(uni->get_quality());
                        duplex_counts[i + n_base_calls] += uni->get_duplex();
                        overlap_counts[i + n_base_calls] += uni->get_overlap();
                        reverse_counts[i + n_base_calls] += uni->get_reverse();
                        qscore_sums[i + n_base_calls] += uni->get_quality();
                    }
                }
                allele_passes[i + n_base_calls] = (duplex_counts[i + n_base_calls] >= aux->conf.min_duplex &&
                                                   counts[i + n_base_calls] >= aux->conf.min_count &&
                                                   overlap_counts[i + n_base_calls] >= aux->conf.min_overlap);
            }
            if((match = tumor.templates.find(base_calls[i])) != tumor.templates.end()) {
                counts[i] = match->second.size();
                for(auto uni: tumor.templates[base_calls[i]]) {
                    if(uni->get_quality() < aux->conf.minPV || uni->get_agreed() < aux->conf.minFA
                            || (float)uni->get_agreed() / uni->get_size() < aux->conf.min_fr) {
                        tsuspect_phreds[i].push_back(uni->get_quality());
                        uni->set_pass(0);
                        ++failed_counts[i];
                    } else {
                        tconfident_phreds[i].push_back(uni->get_quality());
                        duplex_counts[i] += uni->get_duplex();
                        overlap_counts[i] += uni->get_overlap();
                        reverse_counts[i] += uni->get_reverse();
                        qscore_sums[i] += uni->get_quality();
                    }
                }
                allele_passes[i] = (duplex_counts[i] >= aux->conf.min_duplex &&
                                    counts[i] >= aux->conf.min_count &&
                                    overlap_counts[i] >= aux->conf.min_overlap);
                somatic[i] = allele_passes[i] && !allele_passes[i + n_base_calls];
            }
        }
        const int total_depth_tumor(std::accumulate(counts.begin(), counts.begin() + n_base_calls, 0));
        const int total_depth_normal(std::accumulate(counts.begin() + n_base_calls, counts.end(), 0));
        //LOG_DEBUG("Got total depths %i,%i.\n", total_depth_tumor, total_depth_normal);
        std::vector<float> rv_fractions;
        std::vector<float> allele_fractions;
        std::vector<int> quant_est;
        rv_fractions.reserve(reverse_counts.size());
        allele_fractions.reserve(reverse_counts.size());
        quant_est.reserve(reverse_counts.size());
        for(unsigned i = 0; i < n_base_calls; ++i) {
            rv_fractions.push_back((float)counts[i] / reverse_counts[i]);
            allele_fractions.push_back((float)counts[i] / total_depth_tumor);
            quant_est.push_back(estimate_quantity(tconfident_phreds, tsuspect_phreds, i));
        }
        for(unsigned i = 0; i < n_base_calls; ++i) {
            rv_fractions.push_back((float)counts[i + n_base_calls] / reverse_counts[i + n_base_calls]);
            allele_fractions.push_back((float)counts[i + n_base_calls] / total_depth_normal);
            quant_est.push_back(estimate_quantity(nconfident_phreds, nsuspect_phreds, i));
        }
        assert(allele_fractions.size() == 2 * n_base_calls);
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
        bcf_update_format_int32(aux->vcf.vh, vrec, "BMF_QUANT", static_cast<const void *>(quant_est.data()), quant_est.size());
        bcf_update_format_int32(aux->vcf.vh, vrec, "QSS", static_cast<const void *>(qscore_sums.data()), qscore_sums.size());
        bcf_update_format_int32(aux->vcf.vh, vrec, "AMBIG", static_cast<const void *>(ambig), sizeof(ambig));
        bcf_update_format_float(aux->vcf.vh, vrec, "AFR", static_cast<const void *>(allele_fractions.data()), allele_fractions.size());
        bcf_update_info_int32(aux->vcf.vh, vrec, "SOMATIC_CALL", static_cast<const void *>(somatic.data()), somatic.size());
    } /* PairVCFLine::to_bcf */

    void add_hdr_lines(bcf_hdr_t *hdr, const char *lines[], size_t n) {
        while(n--)
            if(bcf_hdr_append(hdr, lines[n]))
                LOG_EXIT("Could not add header line %s. Abort!\n", lines[n]);
    }

    static const char *stack_vcf_lines[] = {
            "##FORMAT=<ID=BMF_PASS,Number=R,Type=Integer,Description=\"1 if variant passes, 0 otherwise.\">",
            "##FORMAT=<ID=ADP,Number=R,Type=Integer,Description=\"Number of unique observations for each allele.\">",
            "##FORMAT=<ID=ADPO,Number=R,Type=Integer,Description=\"Number of unique observations of overlapped read pairs for each allele.\">",
            "##FORMAT=<ID=ADPD,Number=R,Type=Integer,Description=\"Number of duplex observations for each allele. If both reads in an overlapping pair are duplex, this counts each separately.\">",
            "##FORMAT=<ID=ADPR,Number=R,Type=Integer,Description=\"Total number of original reversed reads supporting allele.\">",
            "##FORMAT=<ID=RVF,Number=R,Type=Float,Description=\"Fraction of reads supporting allele which were reversed.\">",
            "##FORMAT=<ID=QSS,Number=R,Type=Integer,Description=\"Q Score Sum for each allele for each sample.\">",
            "##FORMAT=<ID=AMBIG,Number=1,Type=Integer,Description=\"Number of ambiguous (N) base calls at position.\">",
            "##INFO=<ID=SOMATIC_CALL,Number=R,Type=Integer,Description=\"Boolean value for a somatic call for each allele.\">",
            "##INFO=<ID=BMF_QUANT,Number=A,Type=Integer,Description=\"Estimated quantitation for variant allele.\">"
    };

    void add_stack_lines(bcf_hdr_t *hdr) {
        add_hdr_lines(hdr, stack_vcf_lines, COUNT_OF(stack_vcf_lines));
    }


} /* namespace bmf */
