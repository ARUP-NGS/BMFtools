#include "stack.h"

#include <cstdint>
#include <algorithm>
#include "include/igamc_cephes.h"
#include "dlib/misc_util.h"

namespace bmf {
SampleVCFPos::SampleVCFPos(std::unordered_map<std::string, UniqueObservation>& obs, int32_t tid, int32_t pos):
size(obs.size()),
pos(pos),
tid(tid) {
    templates.reserve(4);
    for(auto& pair: obs) templates[pair.second.base_call].push_back(&pair.second);
}

void SampleVCFPos::to_bcf(bcf1_t *vrec, stack_aux_t *aux, const char refbase) {
    int ambig(0);
    unsigned i;
    std::unordered_set<char> base_set{refbase};
    for(auto& pair: templates) {
        if(pair.first == 'N')
           ambig = pair.second.size();
        else
            base_set.insert(pair.first);
    }
    std::vector<char> base_calls(base_set.begin(), base_set.end());
    const size_t n_base_calls(base_calls.size());
    std::vector<std::vector<uint32_t>> confident_phreds;
    std::vector<std::vector<uint32_t>> suspect_phreds;
    confident_phreds.reserve(n_base_calls);
    suspect_phreds.reserve(n_base_calls);
    // Sort lexicographically AFTER putting the reference base first.
    std::sort(base_calls.begin(), base_calls.end(), [refbase](const char a, const char b) {
        return (a == refbase) ? true : (b == refbase) ? false: a < b;
    });
    std::vector<int> counts(n_base_calls);
    std::vector<int> pv_failed(n_base_calls);
    std::vector<int> fr_failed(n_base_calls);
    std::vector<int> fa_failed(n_base_calls);
    std::vector<int> fm_failed(n_base_calls);
    std::vector<int> md_failed(n_base_calls);
    std::vector<int> duplex_counts(n_base_calls);
    std::vector<int> overlap_counts(n_base_calls);
    std::vector<int> reverse_counts(n_base_calls);
    std::vector<int> failed_counts(n_base_calls);
    std::vector<int> allele_passes(n_base_calls);
    std::vector<int> qscore_sums(n_base_calls);
    std::vector<int> rv_sums(n_base_calls);
    vrec->rid = tid;
    vrec->pos = pos;
    vrec->qual = 0;
    vrec->n_sample = 1;
    auto match = templates.find(refbase);
    suspect_phreds.emplace_back();
    confident_phreds.emplace_back();
    if(match != templates.end()) {
        counts[0] = match->second.size();
        for(auto&& uni: match->second) {
            if(uni->get_size() < (unsigned)aux->conf.minFM) {
                uni->pass = 0;
                ++fm_failed[0];
            }
            if(uni->get_quality() < aux->conf.minPV) {
                uni->pass = 0;
                ++pv_failed[0];
            }
            if(uni->get_agreed() < aux->conf.minFA) {
                uni->pass = 0;
                ++fa_failed[0];
            }
            if(aux->conf.md_thresh && uni->md >= aux->conf.md_thresh) {
                uni->pass = 0;
                ++md_failed[0];
            }
            if(uni->get_frac() < aux->conf.min_fr) {
                uni->pass = 0;
                ++fr_failed[0];
            }
            if(!uni->pass) {
                suspect_phreds[0].push_back(uni->get_quality());
                ++failed_counts[0];
            } else {
                confident_phreds[0].push_back(uni->get_quality());
                duplex_counts[0] += uni->get_duplex();
                overlap_counts[0] += uni->get_overlap();
                reverse_counts[0] += uni->get_reverse();
                qscore_sums[0] += uni->get_quality();
                rv_sums[0] += uni->rv;
            }
        }
        allele_passes[0] = (duplex_counts[0] >= aux->conf.min_duplex &&
                            confident_phreds[0].size() >= (unsigned)aux->conf.min_count &&
                            overlap_counts[0] >= aux->conf.min_overlap);
    }

    kstring_t allele_str{0, 0, nullptr};
    ks_resize(&allele_str, 8uL);
    kputc(refbase, &allele_str);
    for(i = 1; i < n_base_calls; ++i) {
        kputc(',', &allele_str), kputc(base_calls[i], &allele_str);
        suspect_phreds.emplace_back();
        confident_phreds.emplace_back();
        if((match = templates.find(base_calls[i])) != templates.end()) {
            counts[i] = match->second.size();
            for(auto&& uni: match->second) {
                if(uni->get_size() < (unsigned)aux->conf.minFM) {
                    uni->pass = 0;
                    ++fm_failed[i];
                }
                if(uni->get_quality() < aux->conf.minPV) {
                    uni->pass = 0;
                    ++pv_failed[i];
                }
                if(uni->get_agreed() < aux->conf.minFA) {
                    uni->pass = 0;
                    ++fa_failed[i];
                }
                if((float)uni->get_agreed() / uni->get_size() < aux->conf.min_fr) {
                    uni->pass = 0;
                    ++fr_failed[i];
                }
                if(aux->conf.md_thresh && uni->md >= aux->conf.md_thresh) {
                    uni->pass = 0;
                    ++md_failed[i];
                }
                if(!uni->pass) {
                    suspect_phreds[i].push_back(uni->get_quality());
                    ++failed_counts[i];
                } else {
                    confident_phreds[i].push_back(uni->get_quality());
                    duplex_counts[i] += uni->get_duplex();
                    overlap_counts[i] += uni->get_overlap();
                    reverse_counts[i] += uni->get_reverse();
                    qscore_sums[i] += uni->get_quality();
                    rv_sums[i] += uni->rv;
                }
            }
            allele_passes[i] = (duplex_counts[i] >= aux->conf.min_duplex &&
                                confident_phreds[i].size() >= (unsigned)aux->conf.min_count &&
                                overlap_counts[i] >= aux->conf.min_overlap);
        }
    }
    assert(n_base_calls == confident_phreds.size());
    assert(n_base_calls == suspect_phreds.size());
    const int total_depth(std::accumulate(counts.begin(), counts.end(), 0));
#if !NDEBUG
    const int total_rv(std::accumulate(reverse_counts.begin(), reverse_counts.end(), 0));
    LOG_DEBUG("total RV: %i.\n", total_rv);
#endif
    std::vector<float> rv_fractions;
    std::vector<float> allele_fractions;
    std::vector<int> quant_est;
    rv_fractions.reserve(n_base_calls);
    allele_fractions.reserve(n_base_calls);
    quant_est.reserve(n_base_calls);
    for(i = 0; i < n_base_calls; ++i) {
        rv_fractions.push_back((float)reverse_counts[i] / counts[i]);
        allele_fractions.push_back((float)counts[i] / total_depth);
        quant_est.push_back(estimate_quantity(confident_phreds, suspect_phreds, i));
    }
    assert(allele_fractions.size() == n_base_calls);
    std::vector<int> adp_pass;
    assert(duplex_counts.size() == n_base_calls);
    adp_pass.reserve(n_base_calls);
    for(auto i: confident_phreds) adp_pass.push_back(static_cast<int>(i.size()));
    bcf_update_alleles_str(aux->vcf.vh, vrec, allele_str.s), free(allele_str.s);
    bcf_int32_vec(aux->vcf.vh, vrec, "ADP_ALL", counts);
    bcf_int32_vec(aux->vcf.vh, vrec, "ADP_PASS", adp_pass);
    bcf_int32_vec(aux->vcf.vh, vrec, "ADPD", duplex_counts);
    bcf_int32_vec(aux->vcf.vh, vrec, "ADPO", overlap_counts);
    bcf_int32_vec(aux->vcf.vh, vrec, "ADPR", reverse_counts);
    bcf_int32_vec(aux->vcf.vh, vrec, "ADPRV", rv_sums);
    bcf_int32_vec(aux->vcf.vh, vrec, "BMF_PASS", allele_passes);
    bcf_int32_vec(aux->vcf.vh, vrec, "BMF_QUANT", quant_est);
    bcf_int32_vec(aux->vcf.vh, vrec, "FA_FAILED", fa_failed);
    bcf_int32_vec(aux->vcf.vh, vrec, "FM_FAILED", fm_failed);
    bcf_int32_vec(aux->vcf.vh, vrec, "FR_FAILED", fr_failed);
    bcf_int32_vec(aux->vcf.vh, vrec, "PV_FAILED", pv_failed);
    bcf_int32_vec(aux->vcf.vh, vrec, "QSS", qscore_sums);
    bcf_update_format_float(aux->vcf.vh, vrec, "REVERSE_FRAC", static_cast<const void *>(rv_fractions.data()), rv_fractions.size());
    bcf_update_format_float(aux->vcf.vh, vrec, "AFR", static_cast<const void *>(allele_fractions.data()), allele_fractions.size());
    bcf_update_format_int32(aux->vcf.vh, vrec, "AMBIG", static_cast<const void *>(&ambig), 1);
    if(aux->conf.md_thresh)
        bcf_update_format_int32(aux->vcf.vh, vrec, "MD_FAILED", static_cast<const void *>(md_failed.data()), md_failed.size());
}

void UniqueObservation::add_obs(const bam_pileup1_t& plp, stack_aux_t *aux) {
    LOG_ASSERT(strcmp(qname.c_str(), bam_get_qname(plp.b)) == 0);
#if !NDEBUG
    for(auto tag: {"PV", "FA"})
        if(!bam_aux_get(plp.b, tag)) LOG_WARNING("Missing tag %s.\n", tag);
#endif
    size += bam_itag(plp.b, "FM");
    base2 = plp_bc(plp);
    cycle2 = dlib::arr_qpos(&plp);
    mq2 = (uint32_t)plp.b->core.qual;
    is_reverse2 = bam_is_rev(plp.b);
    is_overlap = 1;
    rv += (uint32_t)dlib::int_tag_zero(plp.b, "RV");
    if(base2 == base1) {
        discordant = 0;
        agreed += ((uint32_t *)dlib::array_tag(plp.b, "FA"))[cycle2];
        quality = agreed_pvalues(quality, ((uint32_t *)dlib::array_tag(plp.b, "PV"))[cycle2]);
        pvalue = std::pow(10, -0.1 * quality);
    } else if(base1 == 'N') {
        discordant = 0;
        base_call = base2;
        agreed = ((uint32_t *)dlib::array_tag(plp.b, "FA"))[cycle2];
        quality = ((uint32_t *)dlib::array_tag(plp.b, "PV"))[cycle2];
        pvalue = std::pow(10, -0.1 * quality);
    } else if(base2 != 'N') {
        discordant = 1;
        base_call = 'N';
        agreed = 0;
        quality = 0;
        pvalue = 1.;
    }
    const int md2(get_mismatch_density(plp, aux));
    if(md2 > md) md = md2;
}
void PairVCFPos::to_bcf(bcf1_t *vrec, stack_aux_t *aux, int ttid, int tpos) {
    unsigned i;
    const char refbase = aux->get_ref_base(ttid, tpos);
    int ambig[2] = {0, 0};
    std::unordered_set<char> base_set{refbase};
    for(auto&& pair: tumor.templates)
        if(pair.first == 'N')
           ambig[0] = pair.second.size();
        else
            base_set.insert(pair.first);
    for(auto&& pair: normal.templates)
        if(pair.first == 'N')
           ambig[1] = pair.second.size();
        else
            base_set.insert(pair.first);
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
    std::vector<int> md_failed(nbc2);
    std::vector<int> duplex_counts(nbc2);
    std::vector<int> overlap_counts(nbc2);
    std::vector<int> reverse_counts(nbc2);
    std::vector<int> failed_counts(nbc2);
    std::vector<int> allele_passes(nbc2);
    std::vector<int> qscore_sums(nbc2);
    std::vector<int> rv_sums(nbc2);
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
        for(auto&& uni: tumor.templates[refbase]) {
            if(uni->get_size() < (unsigned)aux->conf.minFM) {
                uni->pass = 0;
                ++fm_failed[0];
            }
            if(uni->get_quality() < aux->conf.minPV) {
                uni->pass = 0;
                ++pv_failed[0];
            }
            if(uni->get_agreed() < aux->conf.minFA) {
                uni->pass = 0;
                ++fa_failed[0];
            }
            if(aux->conf.md_thresh && uni->md >= aux->conf.md_thresh) {
                uni->pass = 0;
                ++md_failed[0];
            }
            if(uni->get_frac() < aux->conf.min_fr) {
                uni->pass = 0;
                ++fr_failed[0];
            }
            if(!uni->pass) {
                tsuspect_phreds[0].push_back(uni->get_quality());
                ++failed_counts[0];
            } else {
                tconfident_phreds[0].push_back(uni->get_quality());
                duplex_counts[0] += uni->get_duplex();
                overlap_counts[0] += uni->get_overlap();
                reverse_counts[0] += uni->get_reverse();
                qscore_sums[0] += uni->get_quality();
                rv_sums[0] += uni->rv;
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
        for(auto&& uni: normal.templates[refbase]) {
            if(uni->get_quality() < aux->conf.minPV) {
                uni->pass = 0;
                ++pv_failed[n_base_calls];
            }
            if(uni->get_size() < (unsigned)aux->conf.minFM) {
                uni->pass = 0;
                ++fm_failed[n_base_calls];
            }
            if(aux->conf.md_thresh && uni->md >= aux->conf.md_thresh) {
                uni->pass = 0;
                ++md_failed[n_base_calls];
            }
            if(uni->get_agreed() < aux->conf.minFA) {
                uni->pass = 0;
                ++fa_failed[n_base_calls];
            }
            if((float)uni->get_agreed() / uni->get_size() < aux->conf.min_fr) {
                uni->pass = 0;
                ++fr_failed[n_base_calls];
            }
            if(!uni->pass) {
                nsuspect_phreds[0].push_back(uni->get_quality());
                ++failed_counts[n_base_calls];
            } else {
                nconfident_phreds[0].push_back(uni->get_quality());
                duplex_counts[n_base_calls] += uni->get_duplex();
                overlap_counts[n_base_calls] += uni->get_overlap();
                reverse_counts[n_base_calls] += uni->get_reverse();
                qscore_sums[n_base_calls] += uni->get_quality();
                rv_sums[n_base_calls] += uni->rv;
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
    somatic.push_back(allele_passes[0] & !allele_passes[n_base_calls]);

    kstring_t allele_str{0, 0, nullptr};
    ks_resize(&allele_str, 8uL);
    kputc(refbase, &allele_str);
    for(i = 1; i < n_base_calls; ++i) {
        kputc(',', &allele_str), kputc(base_calls[i], &allele_str);
        nsuspect_phreds.emplace_back();
        nconfident_phreds.emplace_back();
        tsuspect_phreds.emplace_back();
        tconfident_phreds.emplace_back();
        if((match = normal.templates.find(base_calls[i])) != normal.templates.end()) {
            counts[i + n_base_calls] = match->second.size();
            for(auto&& uni: match->second) {
                if(uni->get_size() < (unsigned)aux->conf.minFM) {
                    uni->pass = 0;
                    ++fm_failed[i + n_base_calls];
                }
                if(uni->get_quality() < aux->conf.minPV) {
                    uni->pass = 0;
                    ++pv_failed[i + n_base_calls];
                }
                if(aux->conf.md_thresh && uni->md >= aux->conf.md_thresh) {
                    uni->pass = 0;
                    ++md_failed[i + n_base_calls];
                }
                if(uni->get_agreed() < aux->conf.minFA) {
                    uni->pass = 0;
                    ++fa_failed[i + n_base_calls];
                }
                if((float)uni->get_agreed() / uni->get_size() < aux->conf.min_fr) {
                    uni->pass = 0;
                    ++fr_failed[i + n_base_calls];
                }
                if(!uni->pass) {
                    nsuspect_phreds[i].push_back(uni->get_quality());
                    ++failed_counts[i + n_base_calls];
                } else {
                    nconfident_phreds[i].push_back(uni->get_quality());
                    duplex_counts[i + n_base_calls] += uni->get_duplex();
                    overlap_counts[i + n_base_calls] += uni->get_overlap();
                    reverse_counts[i + n_base_calls] += uni->get_reverse();
                    qscore_sums[i + n_base_calls] += uni->get_quality();
                    rv_sums[i + n_base_calls] += uni->rv;
                }
            }
            allele_passes[i + n_base_calls] = (duplex_counts[i + n_base_calls] >= aux->conf.min_duplex &&
                                               nconfident_phreds[i].size() >= (unsigned)aux->conf.min_count &&
                                               overlap_counts[i + n_base_calls] >= aux->conf.min_overlap);
        }
        if((match = tumor.templates.find(base_calls[i])) != tumor.templates.end()) {
            counts[i] = match->second.size();
            for(auto&& uni: match->second) {
                if(uni->get_quality() < aux->conf.minPV) {
                    uni->pass = 0;
                    ++pv_failed[i];
                }
                if(uni->get_size() < (unsigned)aux->conf.minFM) {
                    uni->pass = 0;
                    ++fm_failed[i];
                }
                if(uni->get_agreed() < aux->conf.minFA) {
                    uni->pass = 0;
                    ++fa_failed[i];
                }
                if((float)uni->get_agreed() / uni->get_size() < aux->conf.min_fr) {
                    uni->pass = 0;
                    ++fr_failed[i];
                }
                if(aux->conf.md_thresh && uni->md >= aux->conf.md_thresh) {
                    uni->pass = 0;
                    ++md_failed[i];
                }
                if(!uni->pass) {
                    tsuspect_phreds[i].push_back(uni->get_quality());
                    ++failed_counts[i];
                } else {
                    tconfident_phreds[i].push_back(uni->get_quality());
                    duplex_counts[i] += uni->get_duplex();
                    overlap_counts[i] += uni->get_overlap();
                    reverse_counts[i] += uni->get_reverse();
                    qscore_sums[i] += uni->get_quality();
                    rv_sums[i] += uni->rv;
                }
            }
            allele_passes[i] = (duplex_counts[i] >= aux->conf.min_duplex &&
                                tconfident_phreds[i].size() >= (unsigned)aux->conf.min_count &&
                                overlap_counts[i] >= aux->conf.min_overlap);
        }
        somatic.push_back(allele_passes[i] & !allele_passes[i + n_base_calls]);
    }
    assert(n_base_calls == tconfident_phreds.size());
    assert(n_base_calls == tsuspect_phreds.size());
    assert(n_base_calls == nconfident_phreds.size());
    assert(n_base_calls == nsuspect_phreds.size());
    const int total_depth_tumor(std::accumulate(counts.begin(), counts.begin() + n_base_calls, 0));
    const int total_depth_normal(std::accumulate(counts.begin() + n_base_calls, counts.end(), 0));
    //LOG_DEBUG("Got total depths %i,%i.\n", total_depth_tumor, total_depth_normal);
    std::vector<float> rv_fractions, allele_fractions;
    std::vector<int> quant_est;
    rv_fractions.reserve(nbc2);
    allele_fractions.reserve(nbc2);
    quant_est.reserve(nbc2);
    for(i = 0; i < n_base_calls; ++i) {
        rv_fractions.push_back((float)reverse_counts[i] / counts[i]);
        allele_fractions.push_back((float)counts[i] / total_depth_tumor);
        quant_est.push_back(estimate_quantity(tconfident_phreds, tsuspect_phreds, i));
    }
    for(i = 0; i < n_base_calls; ++i) {
        rv_fractions.push_back((float)reverse_counts[i + n_base_calls] / counts[i + n_base_calls]);
        allele_fractions.push_back((float)counts[i + n_base_calls] / total_depth_normal);
        quant_est.push_back(estimate_quantity(nconfident_phreds, nsuspect_phreds, i));
    }
    assert(allele_fractions.size() == 2 * n_base_calls);
    std::vector<int> adp_pass;
    assert(duplex_counts.size() == nbc2);
    adp_pass.reserve(nbc2);
    for(auto i: tconfident_phreds) adp_pass.push_back(static_cast<int>(i.size()));
    for(auto i: nconfident_phreds) adp_pass.push_back(static_cast<int>(i.size()));
    bcf_update_alleles_str(aux->vcf.vh, vrec, allele_str.s), free(allele_str.s);
#if !NDEBUG
    if(vrec->pos == 55249070) {
        LOG_DEBUG("Number of base calls: %lu. Size of allele_passes: %lu. Ref: %c\n", n_base_calls, allele_passes.size(), vrec->d.allele[0][0]);
    }
    assert(allele_passes.size() == nbc2);
    assert(fa_failed.size() == nbc2);
#endif
    // @Daniel TODO: Use normal quantity to filter out technical noise.
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
    if(aux->conf.md_thresh)
        bcf_int32_vec(aux->vcf.vh, vrec, "MD_FAILED", md_failed);
    bcf_int32_vec(aux->vcf.vh, vrec, "QSS", qscore_sums);
    bcf_int32_vec(aux->vcf.vh, vrec, "RVC", rv_sums);
    bcf_update_format_float(aux->vcf.vh, vrec, "REVERSE_FRAC", static_cast<const void *>(rv_fractions.data()), rv_fractions.size() * 2);
    bcf_update_format_int32(aux->vcf.vh, vrec, "AMBIG", static_cast<const void *>(ambig), COUNT_OF(ambig) * 2);
    assert(somatic.size() == n_base_calls);
    bcf_update_info_int32(aux->vcf.vh, vrec, "SOMATIC_CALL", static_cast<const void *>(somatic.data()), somatic.size());
} /* PairVCFLine::to_bcf */

static const char *stack_vcf_lines[] = {
        "##INFO=<ID=SOMATIC_CALL,Number=R,Type=Integer,Description=\"Boolean value for a somatic call for each allele.\">",
        "##FORMAT=<ID=ADP_ALL,Number=R,Type=Integer,Description=\"Number of all unique observations for each allele, inc. both low- and high-confidence.\">",
        "##FORMAT=<ID=ADPD,Number=R,Type=Integer,Description=\"Number of duplex observations for each allele. If both reads in an overlapping pair are duplex, this counts each separately.\">",
        "##FORMAT=<ID=ADPO,Number=R,Type=Integer,Description=\"Number of unique observations of overlapped read pairs for each allele.\">",
        "##FORMAT=<ID=ADP_PASS,Number=.,Type=Integer,Description=\"Number of high-confidence unique observations for each allele.\">",
        "##FORMAT=<ID=ADPR,Number=R,Type=Integer,Description=\"Total number of reads aligned to reverse strand.\">",
        "##FORMAT=<ID=ADPRV,Number=R,Type=Integer,Description=\"Number of reads supporting allele which were reversed (inline chemistry).\">",
        "##FORMAT=<ID=AFR,Number=R,Type=Float,Description=\"Allele fractions per allele, including the reference allele.\">",
        "##FORMAT=<ID=AMBIG,Number=1,Type=Integer,Description=\"Number of ambiguous (N) base calls at position.\">",
        "##FORMAT=<ID=BMF_PASS,Number=R,Type=Integer,Description=\"1 if variant passes, 0 otherwise.\">",
        "##FORMAT=<ID=BMF_QUANT,Number=R,Type=Integer,Description=\"Estimated quantitation for each allele.\">",
        "##FORMAT=<ID=AF_FAILED,Number=1,Type=Integer,Description=\"Number of observations failed per sample for aligned fraction below minimm.\">",
        "##FORMAT=<ID=FA_FAILED,Number=R,Type=Integer,Description=\"Number of observations failed per sample for number of supporting observations.\">",
        "##FORMAT=<ID=FM_FAILED,Number=R,Type=Integer,Description=\"Number of observations failed per sample for family size.\">",
        "##FORMAT=<ID=FP_FAILED,Number=1,Type=Integer,Description=\"Number of observations failed per sample for being a barcode QC fail.\">",
        "##FORMAT=<ID=FR_FAILED,Number=R,Type=Integer,Description=\"Number of observations failed per sample for fraction agreed.\">",
        "##FORMAT=<ID=MD_FAILED,Number=R,Type=Integer,Description=\"Number of observations failed for mismatch density filter.\">",
        "##FORMAT=<ID=IMPROPER,Number=1,Type=Integer,Description=\"Number of reads per sample labeled as not being in a proper pair.\">",
        "##FORMAT=<ID=MQ_FAILED,Number=1,Type=Integer,Description=\"Number of observations failed per sample for insufficient mapping quality.\">",
        "##FORMAT=<ID=OVERLAP,Number=1,Type=Integer,Description=\"Number of overlapping read pairs.\">",
        "##FORMAT=<ID=PV_FAILED,Number=R,Type=Integer,Description=\"Number of observations failed per sample for p value cutoff.\">",
        "##FORMAT=<ID=QSS,Number=R,Type=Integer,Description=\"Q Score Sum for each allele for each sample.\">",
        "##FORMAT=<ID=REVERSE_FRAC,Number=R,Type=Float,Description=\"Fraction of reads supporting allele aligned to the reverse strand.\">"
};

void add_stack_lines(bcf_hdr_t *hdr) {
    for(auto line: stack_vcf_lines)
        if(bcf_hdr_append(hdr, line))
            LOG_EXIT("Could not add header line %s. Abort!\n", line);
}


} /* namespace bmf */
