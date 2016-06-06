#include "stack.h"

#include <cstdint>
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
        const size_t nbc2(n_base_calls * 2);
        std::vector<std::vector<uint32_t>> tconfident_phreds;
        std::vector<std::vector<uint32_t>> tsuspect_phreds;
        std::vector<std::vector<uint32_t>> nconfident_phreds;
        std::vector<std::vector<uint32_t>> nsuspect_phreds;
        tconfident_phreds.reserve(n_base_calls);
        tsuspect_phreds.reserve(n_base_calls);
        nconfident_phreds.reserve(n_base_calls);
        nsuspect_phreds.reserve(n_base_calls);
        //TODO: Switch out about two lines for below two.
        /*
        nconfident_phreds.reserve(n_base_calls);
        nsuspect_phreds.reserve(n_base_calls);
        */
        // Sort lexicographically AFTER putting the reference base first.
        std::sort(base_calls.begin(), base_calls.end(), [refbase](const char a, const char b) {
            return (a == refbase) ? true : (b == refbase) ? false: a < b;
        });
        std::vector<int> counts(nbc2);
        std::vector<int> pv_failed(nbc2);
        std::vector<int> fr_failed(nbc2);
        std::vector<int> fa_failed(nbc2);
        std::vector<int> fm_failed(nbc2);
        std::vector<int> duplex_counts(nbc2);
        std::vector<int> overlap_counts(nbc2);
        std::vector<int> reverse_counts(nbc2);
        std::vector<int> failed_counts(nbc2);
        std::vector<int> allele_passes(nbc2);
        std::vector<int> qscore_sums(nbc2);
        std::vector<int> somatic;
        somatic.reserve(n_base_calls);
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
            counts[0] = match->second.size();
            for(auto uni: tumor.templates[refbase]) {
                if(uni->get_size() < (unsigned)aux->conf.minFM) {
                    uni->set_pass(0);
                    ++fm_failed[0];
                }
                if(uni->get_quality() < aux->conf.minPV) {
                    uni->set_pass(0);
                    ++pv_failed[0];
                }
                if(uni->get_agreed() < aux->conf.minFA) {
                    uni->set_pass(0);
                    ++fa_failed[0];
                }
                if(uni->get_frac() < aux->conf.min_fr) {
                    uni->set_pass(0);
                    ++fr_failed[0];
                }
                if(!uni->is_pass()) {
                    tsuspect_phreds[0].push_back(uni->get_quality());
                    ++failed_counts[0];
                } else {
                    tconfident_phreds[0].push_back(uni->get_quality());
                    duplex_counts[0] += uni->get_duplex();
                    overlap_counts[0] += uni->get_overlap();
                    reverse_counts[0] += uni->get_reverse();
                    qscore_sums[0] += uni->get_quality();
                }
            }
            allele_passes[0] = (duplex_counts[0] >= aux->conf.min_duplex &&
                                tconfident_phreds[0].size() >= (unsigned)aux->conf.min_count &&
                                overlap_counts[0] >= aux->conf.min_overlap);
#if !NDEBUG
            if(!allele_passes[0]) {
                LOG_DEBUG("Ref allele failed at %i:%i.\n", vrec->rid, vrec->pos);
                if(duplex_counts[0] < aux->conf.min_duplex) {
                    LOG_DEBUG("Ref allele failed duplex threshold.\n");
                }
                else if(tconfident_phreds[0].size() < (unsigned)aux->conf.min_count) {
                    LOG_DEBUG("Ref allele failed count threshold.\n");
                }
                else if(overlap_counts[0] < aux->conf.min_overlap) {
                    LOG_DEBUG("Failed olap threshold.\n");
                }
            }
#endif
        }
        if((match = normal.templates.find(refbase)) != normal.templates.end()) {
            counts[n_base_calls] = match->second.size();
            for(auto uni: normal.templates[refbase]) {
                if(uni->get_quality() < aux->conf.minPV) {
                    uni->set_pass(0);
                    ++pv_failed[n_base_calls];
                }
                if(uni->get_size() < (unsigned)aux->conf.minFM) {
                    uni->set_pass(0);
                    ++fm_failed[n_base_calls];
                }
                if(uni->get_agreed() < aux->conf.minFA) {
                    uni->set_pass(0);
                    ++fa_failed[n_base_calls];
                }
                if((float)uni->get_agreed() / uni->get_size() < aux->conf.min_fr) {
                    uni->set_pass(0);
                    ++fr_failed[n_base_calls];
                }
                if(!uni->is_pass()) {
                    nsuspect_phreds[0].push_back(uni->get_quality());
                    ++failed_counts[n_base_calls];
                } else {
                    nconfident_phreds[0].push_back(uni->get_quality());
                    duplex_counts[n_base_calls] += uni->get_duplex();
                    overlap_counts[n_base_calls] += uni->get_overlap();
                    reverse_counts[n_base_calls] += uni->get_reverse();
                    qscore_sums[n_base_calls] += uni->get_quality();
                }
            }
            allele_passes[n_base_calls] = (duplex_counts[n_base_calls] >= aux->conf.min_duplex &&
                                           nconfident_phreds[0].size() >= (unsigned)aux->conf.min_count &&
                                           overlap_counts[n_base_calls] >= aux->conf.min_overlap);
#if !NDEBUG
                if(allele_passes[n_base_calls]) {
                    //LOG_DEBUG("Normal reference allele is passing duplex. %i > %i.", duplex_counts[n_base_calls], aux->conf.min_duplex);
                    //LOG_DEBUG("Normal reference allele is passing counts. %lu > %i.", nconfident_phreds[0].size(), aux->conf.min_count);
                    //LOG_DEBUG("Normal allele is overlap counts. %lu > %i.\n", overlap_counts[n_base_calls], aux->conf.min_overlap);
                } else {
                    LOG_DEBUG("Fail!\n Which of the below?\n");
                    LOG_DEBUG("Normal reference allele duplex. %i > %i.", duplex_counts[n_base_calls], aux->conf.min_duplex);
                    LOG_DEBUG("Normal reference allele counts. %lu > %i.", nconfident_phreds[0].size(), aux->conf.min_count);
                    LOG_DEBUG("Normal reference allele overlap counts. %lu > %i.\n", overlap_counts[n_base_calls], aux->conf.min_overlap);
                }
#endif
        }
        somatic.push_back(allele_passes[0] && !allele_passes[n_base_calls]);

        kstring_t allele_str = {0, 0, nullptr};
        ks_resize(&allele_str, 8uL);
        kputc(refbase, &allele_str);
        for(unsigned i = 1; i < n_base_calls; ++i) {
            kputc(',', &allele_str), kputc(base_calls[i], &allele_str);
            nsuspect_phreds.emplace_back();
            nconfident_phreds.emplace_back();
            tsuspect_phreds.emplace_back();
            tconfident_phreds.emplace_back();
            if((match = normal.templates.find(base_calls[i])) != normal.templates.end()) {
                counts[i + n_base_calls] = match->second.size();
                for(auto uni: match->second) {
                    if(uni->get_size() < (unsigned)aux->conf.minFM) {
                        uni->set_pass(0);
                        ++fm_failed[i + n_base_calls];
                    }
                    if(uni->get_quality() < aux->conf.minPV) {
                        uni->set_pass(0);
                        ++pv_failed[i + n_base_calls];
                    }
                    if(uni->get_agreed() < aux->conf.minFA) {
                        uni->set_pass(0);
                        ++fa_failed[i + n_base_calls];
                    }
                    if((float)uni->get_agreed() / uni->get_size() < aux->conf.min_fr) {
                        uni->set_pass(0);
                        ++fr_failed[i + n_base_calls];
                    }
                    if(!uni->is_pass()) {
                        nsuspect_phreds[i].push_back(uni->get_quality());
                        ++failed_counts[i + n_base_calls];
                    } else {
                        nconfident_phreds[i].push_back(uni->get_quality());
                        duplex_counts[i + n_base_calls] += uni->get_duplex();
                        overlap_counts[i + n_base_calls] += uni->get_overlap();
                        reverse_counts[i + n_base_calls] += uni->get_reverse();
                        qscore_sums[i + n_base_calls] += uni->get_quality();
                    }
                }
                allele_passes[i + n_base_calls] = (duplex_counts[i + n_base_calls] >= aux->conf.min_duplex &&
                                                   nconfident_phreds[i].size() >= (unsigned)aux->conf.min_count &&
                                                   overlap_counts[i + n_base_calls] >= aux->conf.min_overlap);
#if 0
                if(allele_passes[i + n_base_calls]) {
                    LOG_DEBUG("Normal allele is passing duplex. %i > %i.", duplex_counts[i + n_base_calls], aux->conf.min_duplex);
                    LOG_DEBUG("Normal allele is passing counts. %lu > %i.", nconfident_phreds[i].size(), aux->conf.min_count);
                    LOG_DEBUG("Normal allele is overlap counts. %lu > %i.\n", overlap_counts[i + n_base_calls], aux->conf.min_overlap);
                } else {
                    LOG_DEBUG("Fail!\n Which of the below?\n");
                    if(duplex_counts[i + n_base_calls] < aux->conf.min_duplex) {
                        LOG_DEBUG("Normal allele duplex. %i > %i.\n", duplex_counts[i + n_base_calls], aux->conf.min_duplex)
                    }
                    if(nconfident_phreds[i].size() < aux->conf.min_count) {
                        LOG_DEBUG("Normal allele counts. %lu > %i. unconfident: %lu.\n", nconfident_phreds[i].size(), aux->conf.min_count,
                                nsuspect_phreds[i].size());
                    }
                    if(overlap_counts[i + n_base_calls] < aux->conf.min_overlap) {
                        LOG_DEBUG("Normal allele overlap counts. %lu > %i.\n", overlap_counts[i + n_base_calls], aux->conf.min_overlap);
                    }
                }
#endif
            }
            if((match = tumor.templates.find(base_calls[i])) != tumor.templates.end()) {
                counts[i] = match->second.size();
                for(auto uni: match->second) {
                    if(uni->get_quality() < aux->conf.minPV) {
                        uni->set_pass(0);
                        ++pv_failed[i];
                    }
                    if(uni->get_size() < (unsigned)aux->conf.minFM) {
                        uni->set_pass(0);
                        ++fm_failed[i];
                    }
                    if(uni->get_agreed() < aux->conf.minFA) {
                        uni->set_pass(0);
                        ++fa_failed[i];
                    }
                    if((float)uni->get_agreed() / uni->get_size() < aux->conf.min_fr) {
                        uni->set_pass(0);
                        ++fr_failed[i];
                    }
                    if(!uni->is_pass()) {
                        tsuspect_phreds[i].push_back(uni->get_quality());
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
                                    tconfident_phreds.size() >= (unsigned)aux->conf.min_count &&
                                    overlap_counts[i] >= aux->conf.min_overlap);
                if(vrec->rid == 6 && vrec->pos == 55249070) {
                    //LOG_DEBUG("counts: %i. min: %i.\n", duplex_counts[i], aux->conf.min_duplex);
                    LOG_DEBUG("%i/%i %s filter for base call %c.\n", vrec->rid, vrec->pos,
                              duplex_counts[i] >= aux->conf.min_duplex ? "pass": "fail", base_calls[i]);
                    LOG_DEBUG("Duplex counts: %i.\n", duplex_counts[i]);
                    //LOG_DEBUG("%i/%i pass duplex filter? %i, (%i, %i).\n", vrec->rid, vrec->pos, (int)(duplex_counts[n_base_calls] >= aux->conf.min_duplex),
                    //          duplex_counts[n_base_calls], aux->conf.min_duplex);
                }
            }
            somatic.push_back(allele_passes[i] && !allele_passes[i + n_base_calls]);
        }
        assert(n_base_calls == tconfident_phreds.size());
        assert(n_base_calls == tsuspect_phreds.size());
        assert(n_base_calls == nconfident_phreds.size());
        assert(n_base_calls == nsuspect_phreds.size());
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
        std::vector<int> adp_pass;
        assert(duplex_counts.size() == nbc2);
        adp_pass.reserve(nbc2);
        for(auto& i: tconfident_phreds) adp_pass.push_back(static_cast<int>(i.size()));
        for(auto& i: nconfident_phreds) adp_pass.push_back(static_cast<int>(i.size()));
        bcf_update_alleles_str(aux->vcf.vh, vrec, allele_str.s), free(allele_str.s);
#if !NDEBUG
        if(vrec->pos == 55249070) {
            LOG_DEBUG("Number of base calls: %lu. Size of allele_passes: %lu.\n", n_base_calls, allele_passes.size());
        }
        LOG_DEBUG("Alleles string: %s.\n", allele_str.s);
        assert(allele_passes.size() == nbc2);
        assert(fa_failed.size() == nbc2);
#endif
        bcf_int32_vec(aux->vcf.vh, vrec, "ADP_ALL", counts);
        bcf_int32_vec(aux->vcf.vh, vrec, "ADP_PASS", adp_pass);
        bcf_int32_vec(aux->vcf.vh, vrec, "ADPD", duplex_counts);
        bcf_int32_vec(aux->vcf.vh, vrec, "ADPO", overlap_counts);
        bcf_int32_vec(aux->vcf.vh, vrec, "ADPR", reverse_counts);
        bcf_update_format_float(aux->vcf.vh, vrec, "AFR", static_cast<const void *>(allele_fractions.data()), allele_fractions.size() * 2);
        bcf_int32_vec(aux->vcf.vh, vrec, "BMF_PASS", allele_passes);
        bcf_int32_vec(aux->vcf.vh, vrec, "BMF_QUANT", quant_est);
        bcf_int32_vec(aux->vcf.vh, vrec, "FA_FAILED", fa_failed);
        bcf_int32_vec(aux->vcf.vh, vrec, "FM_FAILED", fm_failed);
        bcf_int32_vec(aux->vcf.vh, vrec, "FR_FAILED", fr_failed);
        bcf_int32_vec(aux->vcf.vh, vrec, "PV_FAILED", pv_failed);
        bcf_int32_vec(aux->vcf.vh, vrec, "QSS", qscore_sums);
        bcf_update_format_float(aux->vcf.vh, vrec, "RVF", static_cast<const void *>(rv_fractions.data()), rv_fractions.size() * 2);
        bcf_update_format_int32(aux->vcf.vh, vrec, "AMBIG", static_cast<const void *>(ambig), COUNT_OF(ambig) * 2);
        assert(somatic.size() == n_base_calls);
        bcf_update_info_int32(aux->vcf.vh, vrec, "SOMATIC_CALL", static_cast<const void *>(somatic.data()), somatic.size());
    } /* PairVCFLine::to_bcf */

    void add_hdr_lines(bcf_hdr_t *hdr, const char *lines[], size_t n) {
        while(n--)
            if(bcf_hdr_append(hdr, lines[n]))
                LOG_EXIT("Could not add header line %s. Abort!\n", lines[n]);
    }

    static const char *stack_vcf_lines[] = {
            "##INFO=<ID=SOMATIC_CALL,Number=R,Type=Integer,Description=\"Boolean value for a somatic call for each allele.\">",
            "##FORMAT=<ID=ADP_ALL,Number=R,Type=Integer,Description=\"Number of all unique observations for each allele, inc. both low- and high-confidence.\">",
            "##FORMAT=<ID=ADPD,Number=R,Type=Integer,Description=\"Number of duplex observations for each allele. If both reads in an overlapping pair are duplex, this counts each separately.\">",
            "##FORMAT=<ID=ADPO,Number=R,Type=Integer,Description=\"Number of unique observations of overlapped read pairs for each allele.\">",
            "##FORMAT=<ID=ADP_PASS,Number=.,Type=Integer,Description=\"Number of high-confidence unique observations for each allele.\">",
            "##FORMAT=<ID=ADPR,Number=R,Type=Integer,Description=\"Total number of original reversed reads supporting allele.\">",
            "##FORMAT=<ID=AFR,Number=R,Type=Float,Description=\"Allele fractions per allele, including the reference allele.\">",
            "##FORMAT=<ID=AMBIG,Number=1,Type=Integer,Description=\"Number of ambiguous (N) base calls at position.\">",
            "##FORMAT=<ID=BMF_PASS,Number=R,Type=Integer,Description=\"1 if variant passes, 0 otherwise.\">",
            "##FORMAT=<ID=BMF_QUANT,Number=.,Type=Integer,Description=\"Estimated quantitation for variant allele.\">",
            "##FORMAT=<ID=AF_FAILED,Number=1,Type=Integer,Description=\"Number of observations failed per sample for aligned fraction below minimm.\">",
            "##FORMAT=<ID=FA_FAILED,Number=1,Type=Integer,Description=\"Number of observations failed per sample for number of supporting observations.\">",
            "##FORMAT=<ID=FM_FAILED,Number=1,Type=Integer,Description=\"Number of observations failed per sample for family size.\">",
            "##FORMAT=<ID=FP_FAILED,Number=1,Type=Integer,Description=\"Number of observations failed per sample for being a barcode QC fail.\">",
            "##FORMAT=<ID=FR_FAILED,Number=1,Type=Integer,Description=\"Number of observations failed per sample for fraction agreed.\">",
            "##FORMAT=<ID=IMPROPER,Number=1,Type=Integer,Description=\"Number of reads per sample labeled as not being in a proper pair.\">",
            "##FORMAT=<ID=MQ_FAILED,Number=1,Type=Integer,Description=\"Number of observations failed per sample for insufficient mapping quality.\">",
            "##FORMAT=<ID=OVERLAP,Number=1,Type=Integer,Description=\"Number of overlapping read pairs.\">",
            "##FORMAT=<ID=PV_FAILED,Number=1,Type=Integer,Description=\"Number of observations failed per sample for p value cutoff.\">",
            "##FORMAT=<ID=QSS,Number=R,Type=Integer,Description=\"Q Score Sum for each allele for each sample.\">",
            "##FORMAT=<ID=RVF,Number=R,Type=Float,Description=\"Fraction of reads supporting allele which were reversed.\">",
    };

    void add_stack_lines(bcf_hdr_t *hdr) {
        add_hdr_lines(hdr, stack_vcf_lines, COUNT_OF(stack_vcf_lines));
    }


} /* namespace bmf */
