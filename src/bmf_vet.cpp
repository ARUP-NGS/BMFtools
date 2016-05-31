#include "bmf_vet.h"

#include <getopt.h>
#include <algorithm>
#include <assert.h>

#include "dlib/bam_util.h"
#include "dlib/vcf_util.h"
#include "include/igamc_cephes.h"
#include "htslib/tbx.h"

namespace bmf {

    KHASH_MAP_INIT_STR(names, const bam_pileup1_t *)
    static int max_depth = (1 << 20); // 262144

    void vetter_usage(int retcode)
    {
        fprintf(stderr,
                        "Checks variants from the bcf within a set of bed regions."
                        "bmftools vet does so by examining the bmftools metadata present in the bam.\n"
                        "Usage:\nbmftools vet -o <out.vcf [stdout]> <in.bcf> <in.srt.indexed.bam>\n"
                        "Optional arguments:\n"
                        "-b, --bedpath\tPath to bed file to only validate variants in said region. REQUIRED.\n"
                        "-c, --min-count\tMinimum number of observations for a given allele passing filters to pass variant. Default: 1.\n"
                        "-s, --min-family-size\tMinimum number of reads in a family to include a that collapsed observation\n"
                        "-D, --min-duplex\tMinimum number of duplex reads required to pass a variant.\n"
                        "-2, --skip-secondary\tSkip secondary alignments.\n"
                        "-S, --skip-supplementary\tSkip supplementary alignments.\n"
                        "-q, --skip-qcfail\tSkip reads marked as QC fail.\n"
                        "-r, --skip-duplicates\tSkip reads marked as duplicates.\n"
                        "-P, --skip-improper\tSkip reads not marked as being in a proper pair.\n"
                        "-F, --skip-recommended.\tSkip secondary, qcfail, and pcr duplicates.\n"
                        "-d, --max-depth\tMaximum depth for pileups. Default: %i.\n"
                        "-f, --min-fraction-agreed\tMinimum fraction of reads in a family agreed on a base call\n"
                        "-v, --min-phred-quality\tMinimum calculated p-value on a base call in phred space\n"
                        "-p, --padding\tNumber of bases outside of bed region to pad.\n"
                        "-a, --min-family-agreed\tMinimum number of reads in a family agreed on a base call\n"
                        "-m, --min-mapping-quality\tMinimum mapping quality for reads for inclusion\n"
                        "-B, --emit-bcf-format\tEmit bcf-formatted output. (Defaults to vcf).\n",
                max_depth
                );
        exit(retcode);
    }

    struct vetter_aux_t {
        samFile *fp;
        hts_itr_t *iter;
        bam_hdr_t *header;
        vcfFile *vcf_fp;
        vcfFile *vcf_ofp;
        bcf_hdr_t *vcf_header;
        khash_t(bed) *bed;
        float min_fr; // Minimum fraction of family members agreed on base
        float minAF; // Minimum aligned fraction
        int max_depth;
        int minFM;
        uint32_t minFA;
        uint32_t minPV;
        int min_count;
        int min_duplex;
        int min_overlap;
        uint32_t skip_improper:1;
        uint32_t vet_all:1;
        uint32_t minmq:8;
        uint32_t skip_flag; // Skip reads with any bits set to true
    };


    int vet_core_nobed(vetter_aux_t *aux);
    void vetter_error(const char *message, int retcode)
    {
        fprintf(stderr, message);
        exit(retcode);
    }


    static int read_bam(void *data, bam1_t *b)
    {
        vetter_aux_t *aux = (vetter_aux_t*)data; // data in fact is a pointer to an auxiliary structure
        int ret;
        for(;;)
        {
            if(!aux->iter) LOG_EXIT("Need to access bam with index.\n");
            ret = sam_itr_next(aux->fp, aux->iter, b);
            if ( ret<0 ) break;
            // Skip unmapped, secondary, qcfail, duplicates.
            // Skip improper if option set
            // Skip MQ < minmq
            // Skip FM < minFM
            // Skip AF < minAF
            if ((b->core.flag & aux->skip_flag) ||
                (aux->skip_improper && ((b->core.flag & BAM_FPROPER_PAIR) == 0)) || // Skip improper if set.
                (int)b->core.qual < aux->minmq ||
                (bam_itag(b, "FP") == 0) || (aux->minAF && bam_aux2f(bam_aux_get(b, "AF")) < aux->minAF))
                    continue;
            break;
        }
        return ret;
    }


    /*
     * :param: [bcf1_t *] vrec - Variant record to test.
     *   # UniObs passing
     * 2. What do I want for INFO?
     *
     * Add 'fm' tag to note which families of reads have already had their fm adjusted.
     * Separate from the upper-case tag!
     */
    void bmf_var_tests(bcf1_t *vrec, const bam_pileup1_t *plp, int n_plp, vetter_aux_t *aux, std::vector<int>& pass_values,
            std::vector<int>& n_obs, std::vector<int>& n_duplex, std::vector<int>& n_overlaps, std::vector<int> &n_failed,
            std::vector<int>& quant_est, int& n_all_overlaps, int& n_all_duplex, int& n_all_disagreed) {
        int khr, s, s2, i;
        n_all_disagreed = n_all_overlaps = 0;
        khiter_t k;
        uint32_t *FA1, *PV1, *FA2, *PV2;
        char *qname;
        uint8_t *seq, *seq2, *tmptag;
        std::vector<std::vector<uint32_t>> confident_phreds;
        std::vector<std::vector<uint32_t>> suspect_phreds;
        confident_phreds.reserve(vrec->n_allele);
        suspect_phreds.reserve(vrec->n_allele);
        // Build overlap hash
        khash_t(names) *hash(kh_init(names));
        const int sk = 1;
        // Set the r1/r2 flags for the reads to ignore to 0
        // Set the ones where we see it twice to (BAM_FREAD1 | BAM_FREAD2).
        for(i = 0; i < n_plp; ++i) {
            if(plp[i].is_del || plp[i].is_refskip) continue;
            // Skip any reads failed for FA < minFA or FR < min_fr
            qname = bam_get_qname(plp[i].b);
            k = kh_get(names, hash, qname);
            if(k == kh_end(hash)) {
                k = kh_put(names, hash, qname, &khr);
                kh_val(hash, k) = &plp[i];
            } else {
                ++n_all_overlaps;
                bam_aux_append(plp[i].b, "SK", 'i', sizeof(int), (uint8_t *)&sk); // Skip
                bam_aux_append(kh_val(hash, k)->b, "KR", 'i', sizeof(int), (uint8_t *)&sk); // Keep Read
                if((tmptag = bam_aux_get(kh_val(hash, k)->b, "fm")) == nullptr) {
                    uint8_t *FM1 = bam_aux_get(kh_val(hash, k)->b, "FM");
                    const int FM_sum = bam_aux2i(FM1) + bam_itag(plp[i].b, "FM");
                    bam_aux_del(kh_val(hash, k)->b, FM1);
                    bam_aux_append(kh_val(hash, k)->b, "FM", 'i', sizeof(int), (uint8_t *)&FM_sum);
                    bam_aux_append(kh_val(hash, k)->b, "fm", 'i', sizeof(int), (uint8_t *)&sk);
                    bam_aux_append(plp[i].b, "fm", 'i', sizeof(int), (uint8_t *)&sk);
                }
                PV1 = (uint32_t *)dlib::array_tag(kh_val(hash, k)->b, "PV");
                FA1 = (uint32_t *)dlib::array_tag(kh_val(hash, k)->b, "FA");
                seq = bam_get_seq(kh_val(hash, k)->b);
                s = bam_seqi(seq, kh_val(hash, k)->qpos);
                PV2 = (uint32_t *)dlib::array_tag(plp[i].b, "PV");
                FA2 = (uint32_t *)dlib::array_tag(plp[i].b, "FA");
                seq2 = bam_get_seq(plp[i].b);
                s2 = bam_seqi(seq2, plp[i].qpos);
                const int32_t arr_qpos1 = dlib::arr_qpos(kh_val(hash, k));
                const int32_t arr_qpos2 = dlib::arr_qpos(&plp[i]);
                if(s == s2) {
                    PV1[arr_qpos1] = agreed_pvalues(PV1[arr_qpos1], PV2[arr_qpos2]);
                    FA1[arr_qpos1] = FA1[arr_qpos1] + FA2[arr_qpos2];
                } else if(s == dlib::htseq::HTS_N) {
                    set_base(seq, seq_nt16_str[bam_seqi(seq2, plp[i].qpos)], kh_val(hash, k)->qpos);
                    PV1[arr_qpos1] = PV2[arr_qpos2];
                    FA1[arr_qpos1] = FA2[arr_qpos2];
                } else if(s2 != dlib::htseq::HTS_N) {
                    ++n_all_disagreed;
                    // Disagreed, both aren't N: N the base, set agrees and p values to 0!
                    n_base(seq, kh_val(hash, k)->qpos); // if s2 == dlib::htseq::HTS_N, do nothing.
                    PV1[arr_qpos1] = 0u;
                    FA1[arr_qpos1] = 0u;
                }
            }
        }
        // Reads in the pair have now been merged, and those to be skipped have been tagged "SK".
        for(unsigned j = 0; j < vrec->n_allele; ++j) {
            confident_phreds.emplace_back(); // Make the vector for PVs for this allele.
            suspect_phreds.emplace_back(); // Make the vector for PVs for this allele.
            if(strcmp(vrec->d.allele[j], "<*>") == 0) {
                LOG_DEBUG("Allele is meaningless/useless <*>. Continuing.\n");
                continue;
            }
            const char allele(vrec->d.allele[j][0]);
            //LOG_DEBUG("Checking stack for reads with base %c.\n", allele);
            for(int i = 0; i < n_plp; ++i) {
                if(plp[i].is_del || plp[i].is_refskip) continue;
                if((tmptag = bam_aux_get(plp[i].b, "SK")) != nullptr) continue;

                seq = bam_get_seq(plp[i].b);
                FA1 = (uint32_t *)dlib::array_tag(plp[i].b, "FA");
                PV1 = (uint32_t *)dlib::array_tag(plp[i].b, "PV");
                if(bam_seqi(seq, plp[i].qpos) == seq_nt16_table[(uint8_t)allele]) { // Match!
                    //LOG_DEBUG("Found read supporting allele '%i', '%c'.\n", bam_seqi(seq, plp[i].qpos), allele);
                    const int32_t arr_qpos1 = dlib::arr_qpos(&plp[i]);
                    if(bam_itag(plp[i].b, "FM") < aux->minFM ||
                       FA1[arr_qpos1] < aux->minFA ||
                            PV1[arr_qpos1] < aux->minPV ||
                            (double)FA1[arr_qpos1] / (tmptag ? bam_aux2i(tmptag): 1 ) < aux->min_fr) {
                        /*
                        LOG_DEBUG("Failed a read because FM ? %s FA? %s PV? %s\n",
                                (tmptag ? bam_aux2i(tmptag): 1) < aux->minFM ? "true": "false",
                                FA1[arr_qpos1] < aux->minFA ? "true": "false",
                                PV1[arr_qpos1] < aux->minPV ? "true": "false"
                                );
                        */
                        ++n_failed[j];
                        suspect_phreds[j].push_back(PV1[arr_qpos1]);
                    } else {
                        //LOG_DEBUG("Passed a read!\n");
                        // TODO: replace n_obs vector with calls to size
                        confident_phreds[j].push_back(PV1[arr_qpos1]);
                        ++n_obs[j];
                        if((tmptag = bam_aux_get(plp[i].b, "DR")) != nullptr) {
                            if(bam_aux2i(tmptag)) {
                                ++n_duplex[j]; // Has DR tag and its value is nonzero.
                                //LOG_DEBUG("Found a duplex read!\n");
                            }
                        }
                        if((tmptag = bam_aux_get(plp[i].b, "KR")) != nullptr) {
                            ++n_overlaps[j];
                            bam_aux_del(plp[i].b, tmptag);
                        }
                    }
                }
            }
            pass_values[j] = n_obs[j] >= aux->min_count && n_duplex[j] >= aux->min_duplex && n_overlaps[j] >= aux->min_overlap;
            quant_est[j] = estimate_quantity(confident_phreds, suspect_phreds, j);
            //LOG_DEBUG("Allele #%i pass? %s\n", j + 1, pass_values[j] ? "True": "False");
        }
        // Now estimate the fraction likely correct.
        for(i = 0; i < n_plp; ++i)
            if((tmptag = bam_aux_get(plp[i].b, "SK")) != nullptr)
                bam_aux_del(plp[i].b, tmptag);
        kh_destroy(names, hash);
        n_all_duplex = std::accumulate(n_duplex.begin(), n_duplex.begin() + vrec->n_allele, 0);
    }

    int read_bcf(vetter_aux_t *aux, hts_itr_t *vcf_iter, bcf1_t *vrec)
    {
        return vcf_iter ? bcf_itr_next(aux->vcf_fp, vcf_iter, vrec)
                        : bcf_read1(aux->vcf_fp, aux->vcf_header, vrec);
    }

    int vet_core_bed(vetter_aux_t *aux) {
        int n_plp;
        const bam_pileup1_t *plp(nullptr);
        tbx_t *vcf_idx(nullptr);
        hts_idx_t *bcf_idx(nullptr);
        hts_idx_t *idx(sam_index_load(aux->fp, aux->fp->fn));
        switch(hts_get_format(aux->vcf_fp)->format) {
        case vcf:
            return vet_core_nobed(aux);
            // Tabix indexed vcfs don't work currently. Throws a segfault in htslib.
        case bcf:
            bcf_idx = bcf_index_load(aux->vcf_fp->fn);
            if(!bcf_idx) LOG_EXIT("Could not load CSI index: %s\n", aux->vcf_fp->fn);
            break;
        default:
            LOG_EXIT("Unrecognized variant file type! (%i).\n", hts_get_format(aux->vcf_fp)->format);
            break; // This never happens -- LOG_EXIT exits.
        }
        /*
        if(!(vcf_idx || bcf_idx)) {
            LOG_EXIT("Require an indexed variant file. Abort!\n");
        }
        */
        bcf1_t *vrec = bcf_init();
        // Unpack all shared data -- up through INFO, but not including FORMAT
        vrec->max_unpack = BCF_UN_FMT;
        vrec->rid = -1;
        hts_itr_t *vcf_iter = nullptr;

        std::vector<int32_t> pass_values(NUM_PREALLOCATED_ALLELES);
        std::vector<int32_t> uniobs_values(NUM_PREALLOCATED_ALLELES);
        std::vector<int32_t> duplex_values(NUM_PREALLOCATED_ALLELES);
        std::vector<int32_t> overlap_values(NUM_PREALLOCATED_ALLELES);
        std::vector<int32_t> fail_values(NUM_PREALLOCATED_ALLELES);
        std::vector<int32_t> quant_est(NUM_PREALLOCATED_ALLELES);
        std::vector<khiter_t> keys(dlib::make_sorted_keys(aux->bed));
        for(khiter_t ki: keys) {
            for(unsigned j = 0; j < kh_val(aux->bed, ki).n; ++j) {
                int tid, start, stop, pos = -1;

                // Handle coordinates
                tid = kh_key(aux->bed, ki);
                // rid is set to -1 before use. This won't be triggered.
                start = get_start(kh_val(aux->bed, ki).intervals[j]);
                stop = get_stop(kh_val(aux->bed, ki).intervals[j]);
                //LOG_DEBUG("Beginning to work through region #%i on contig %s:%i-%i.\n", j + 1, aux->header->target_name[tid], start, stop);

                // Fill vcf_iter from tbi or csi index. If both are null, go through the full file.
                vcf_iter = vcf_idx ? tbx_itr_queryi(vcf_idx, tid, start, stop)
                                   :bcf_idx ? bcf_itr_queryi(bcf_idx, tid, start, stop)
                                            : nullptr;
                //vcf_iter = vcf_idx ? hts_itr_query(vcf_idx->idx, tid, start, stop, tbx_readrec): bcf_idx ? bcf_itr_queryi(bcf_idx, tid, start, stop): nullptr;

                int n_disagreed = 0;
                int n_overlapped = 0;
                int n_duplex = 0;
                bam_plp_t pileup = bam_plp_init(read_bam, (void *)aux);
                bam_plp_set_maxcnt(pileup, max_depth);
                if (aux->iter) hts_itr_destroy(aux->iter);
                aux->iter = sam_itr_queryi(idx, tid, start - 500, stop);
                while(read_bcf(aux, vcf_iter, vrec) >= 0) {
                    //if (aux->iter) hts_itr_destroy(aux->iter);
                    if(!bcf_is_snp(vrec)) {
                        LOG_DEBUG("Variant isn't a snp. Skip!\n");
                        bcf_write(aux->vcf_ofp, aux->vcf_header, vrec);
                        continue; // Only handle simple SNVs
                    }
                    if(!dlib::vcf_bed_test(vrec, aux->bed) && !aux->vet_all) {
                        LOG_DEBUG("Outside of bed region. Skip.\n");
                        continue; // Only handle simple SNVs
                    }
                    /*
                    while(j < kh_val(aux->bed, ki).n && vrec->pos > get_stop(kh_val(aux->bed, ki).intervals[j])) {
                        LOG_DEBUG("New interval!\n");
                        ++j;
                    }
                    while(kh_key(aux->bed, ki) < vrec->rid) {
                        LOG_DEBUG("New interval!\n");
                        ki = keys[++kint];
                    }
                    while(j < kh_val(aux->bed, ki).n && vrec->pos > get_stop(kh_val(aux->bed, ki).intervals[j])) ++j;
                    */
                    bcf_unpack(vrec, BCF_UN_STR); // Unpack the allele fields
                    if(pos > vrec->pos) LOG_EXIT("pos is after variant. WTF?");
                    while(pos < vrec->pos && tid <= vrec->rid && (plp = bam_plp_auto(pileup, &tid, &pos, &n_plp)) > 0);
                    //assert(pos == vrec->pos);
                        //LOG_DEBUG("Now at position %i on contig %s with %i pileups.\n", pos, aux->header->target_name[tid], n_plp);
                        /* Zoom ahead until you're at the correct position */
                    if(!plp) {
                        if(n_plp == -1) {
                            LOG_WARNING("Could not make pileup for region %s:%i-%i. n_plp: %i, pos%i, tid%i.\n",
                                        aux->header->target_name[tid], start, stop, n_plp, pos, tid);
                        }
                        else if(n_plp == 0){
                            LOG_WARNING("No reads at position. Skip this variant.\n");
                        } else LOG_EXIT("No pileup stack, but n_plp doesn't signal an error or an empty stack?\n");
                    }
                    //LOG_DEBUG("tid: %i. rid: %i. var pos: %i.\n", tid, vrec->rid, vrec->pos);
                    if(pos != vrec->pos || tid != vrec->rid) {
                        //LOG_DEBUG("BAM: pos: %i. Contig: %s.\n", pos, aux->header->target_name[tid]);
                        LOG_WARNING("Position %s:%i (1-based) not found in pileups in bam. Writing unmodified. Super weird...\n",
                                    aux->header->target_name[vrec->rid], vrec->pos + 1);
                        bcf_write(aux->vcf_ofp, aux->vcf_header, vrec);
                        continue;
                    }
                    // Reset vectors for each pass.
                    memset(uniobs_values.data(), 0, sizeof(int32_t) * uniobs_values.size());
                    memset(duplex_values.data(), 0, sizeof(int32_t) * duplex_values.size());
                    memset(fail_values.data(), 0, sizeof(int32_t) * fail_values.size());
                    memset(overlap_values.data(), 0, sizeof(int32_t) * overlap_values.size());
                    // Perform tests to provide the results for the tags.
                    bmf_var_tests(vrec, plp, n_plp, aux, pass_values, uniobs_values, duplex_values, overlap_values,
                                 fail_values, quant_est, n_overlapped, n_duplex, n_disagreed);
                    // Add tags
                    bcf_update_info_int32(aux->vcf_header, vrec, "DISC_OVERLAP", (void *)&n_disagreed, 1);
                    bcf_update_info_int32(aux->vcf_header, vrec, "OVERLAP", (void *)&n_overlapped, 1);
                    bcf_update_info_int32(aux->vcf_header, vrec, "DUPLEX_DEPTH", (void *)&n_duplex, 1);
                    bcf_update_info(aux->vcf_header, vrec, "BMF_VET", (void *)(&pass_values[0]), vrec->n_allele, BCF_HT_INT);
                    bcf_update_info(aux->vcf_header, vrec, "BMF_FAIL", (void *)(&fail_values[0]), vrec->n_allele, BCF_HT_INT);
                    bcf_update_info(aux->vcf_header, vrec, "BMF_DUPLEX", (void *)&duplex_values[0], vrec->n_allele, BCF_HT_INT);
                    bcf_update_info(aux->vcf_header, vrec, "BMF_UNIOBS", (void *)&uniobs_values[0], vrec->n_allele, BCF_HT_INT);

                    // Pass or fail them individually.
                    bcf_write(aux->vcf_ofp, aux->vcf_header, vrec);
                    //bam_plp_reset(pileup);
                }
                if(vcf_iter) hts_itr_destroy(vcf_iter);
                bam_plp_destroy(pileup);
            }
        }
        if(bcf_idx) hts_idx_destroy(bcf_idx);
        if(vcf_idx) tbx_destroy(vcf_idx);
        if(aux->iter) hts_itr_destroy(aux->iter);
        hts_idx_destroy(idx);
        bcf_destroy(vrec);
        return EXIT_SUCCESS;
    }

    int vet_core_nobed(vetter_aux_t *aux) {
#if !NDEBUG
        int n_skipped = 0;
#endif
        int n_plp;
        const bam_pileup1_t *plp;
        tbx_t *vcf_idx = nullptr;
        hts_idx_t *bcf_idx = nullptr;
        hts_idx_t *idx = sam_index_load(aux->fp, aux->fp->fn);
        switch(hts_get_format(aux->vcf_fp)->format) {
        case vcf:
            //LOG_WARNING("Somehow, tabix reading doesn't seem to work. I'm deleting this index and iterating through the whole vcf.\n");
            vcf_idx = nullptr;
            /*
            if(vcf_idx) {
                tbx_destroy(vcf_idx);
                vcf_idx = nullptr;
            }
            */
            break;
        case bcf:
            bcf_idx = nullptr;
            break;
        default:
            LOG_EXIT("Unrecognized variant file type! (%i).\n", hts_get_format(aux->vcf_fp)->format);
            break; // This never happens -- LOG_EXIT exits.
        }
        bcf1_t *vrec = bcf_init();
        // Unpack all shared data -- up through INFO, but not including FORMAT
        vrec->max_unpack = BCF_UN_FMT;
        vrec->rid = -1;
        hts_itr_t *vcf_iter = nullptr;

        std::vector<int32_t> pass_values(NUM_PREALLOCATED_ALLELES);
        std::vector<int32_t> uniobs_values(NUM_PREALLOCATED_ALLELES);
        std::vector<int32_t> duplex_values(NUM_PREALLOCATED_ALLELES);
        std::vector<int32_t> overlap_values(NUM_PREALLOCATED_ALLELES);
        std::vector<int32_t> fail_values(NUM_PREALLOCATED_ALLELES);
        std::vector<int32_t> quant_est(NUM_PREALLOCATED_ALLELES);
        bam_plp_t pileup(nullptr);

        for(int i = 0; i < aux->header->n_targets; ++i) {

            int pos = -1;
            int tid = i;
            const int start = 0;
            const int stop = aux->header->target_len[tid];

            // Handle coordinates
            // rid is set to -1 before use. This won't be triggered.
            //LOG_DEBUG("Beginning to work through region #%i on contig %s:%i-%i.\n", tid + 1, aux->header->target_name[tid], start, stop);

            // Fill vcf_iter from tbi or csi index. If both are null, go through the full file.

            int n_disagreed = 0;
            int n_overlapped = 0;
            int n_duplex = 0;
            while(read_bcf(aux, vcf_iter, vrec) >= 0) {
                if(!vcf_iter && vrec->rid != tid) {
                    i = tid = vrec->rid; // Finished the last contig, on the next one.
                }
                if(!bcf_is_snp(vrec)) {
                    LOG_DEBUG("Variant isn't a snp. Skip!\n");
                    bcf_write(aux->vcf_ofp, aux->vcf_header, vrec);
                    continue; // Only handle simple SNVs
                }
                bcf_unpack(vrec, BCF_UN_STR); // Unpack the allele fields
                LOG_DEBUG("Querying for tid and pos %i, %i.\n", vrec->rid, vrec->pos);
                aux->iter = sam_itr_queryi(idx, vrec->rid, vrec->pos - 500, vrec->pos);
                //LOG_DEBUG("Before plp_auto tid %i and pos %i for a variant at %i, %i\n", tid, pos, vrec->rid, vrec->pos);
                if(pileup) bam_plp_destroy(pileup);
                pileup = bam_plp_init(read_bam, (void *)aux), bam_plp_set_maxcnt(pileup, max_depth);
                LOG_DEBUG("Max depth: %i.\n", max_depth);
                bam_plp_reset(pileup);
                plp = bam_plp_auto(pileup, &tid, &pos, &n_plp);
                //LOG_DEBUG("Hey, I'm evaluating a variant record now with tid %i and pos %i for a variant at %i, %i\n", tid, pos, vrec->rid, vrec->pos);
                assert(tid == vrec->rid);
                while (pos < vrec->pos &&
                       ((plp = bam_plp_auto(pileup, &tid, &pos, &n_plp)) > 0)) {
                    assert(tid == vrec->rid);
                    LOG_DEBUG("Now at position %i on contig %s with %i pileups.\n", pos, aux->header->target_name[tid], n_plp);
                    /* Zoom ahead until you're at the correct position */
                }
                assert(tid == vrec->rid);
                if(pos != vrec->pos) {
                    LOG_DEBUG("Try again\n");
                    bam_plp_reset(pileup);
                    while (pos < vrec->pos &&
                           ((plp = bam_plp_auto(pileup, &tid, &pos, &n_plp)) > 0)) {
                        assert(tid == vrec->rid);
                        LOG_DEBUG("Now at position %i on contig %s with %i pileups.\n", pos, aux->header->target_name[tid], n_plp);
                        /* Zoom ahead until you're at the correct position */
                    }
                    LOG_EXIT("Pos is not correct (%i rather than expected %i). n_plp: %i.\n", pos, vrec->pos, n_plp);
                }
                //LOG_DEBUG("Hey, I'm evaluating a variant record now with tid %i and pos %i for a variant at %i, %i\n", tid, pos, vrec->rid, vrec->pos);
                if(!plp) {
                    if(n_plp == -1) {
                        LOG_WARNING("Could not make pileup for region %s:%i-%i. n_plp: %i, pos%i, tid%i.\n",
                                    aux->header->target_name[tid], start, stop, n_plp, pos, tid);
                    }
                    else if(n_plp == 0) {
                        //LOG_DEBUG("Vrec rid: %i. Vrec pos: %i. tid: %i. pos: %i.\n", vrec->rid, vrec->pos, tid, pos);
                        LOG_WARNING("No reads at position. Skip this variant.\n");
#if !NDEBUG
                        ++n_skipped;
#endif
                        continue;
                    } else LOG_EXIT("No pileup stack, but n_plp doesn't signal an error or an empty stack?\n");
                }
                //LOG_DEBUG("tid: %i. rid: %i. var pos: %i.\n", tid, vrec->rid, vrec->pos);
                if(pos != vrec->pos || tid != vrec->rid) {
                    //LOG_DEBUG("BAM: pos: %i. Contig: %s.\n", pos, aux->header->target_name[tid]);
                    LOG_WARNING("Position %s:%i (1-based) not found in pileups in bam. Writing unmodified. Super weird...\n", aux->header->target_name[vrec->rid], vrec->pos + 1);
                    bcf_write(aux->vcf_ofp, aux->vcf_header, vrec);
                    continue;
                }
                // Reset vectors for each pass.
                memset(uniobs_values.data(), 0, sizeof(int32_t) * uniobs_values.size());
                memset(duplex_values.data(), 0, sizeof(int32_t) * duplex_values.size());
                memset(fail_values.data(), 0, sizeof(int32_t) * fail_values.size());
                memset(overlap_values.data(), 0, sizeof(int32_t) * overlap_values.size());
                // Perform tests to provide the results for the tags.
                bmf_var_tests(vrec, plp, n_plp, aux, pass_values, uniobs_values, duplex_values, overlap_values,
                             fail_values, quant_est, n_overlapped, n_duplex, n_disagreed);
                // Add tags
                bcf_update_info_int32(aux->vcf_header, vrec, "DISC_OVERLAP", (void *)&n_disagreed, 1);
                bcf_update_info_int32(aux->vcf_header, vrec, "OVERLAP", (void *)&n_overlapped, 1);
                bcf_update_info_int32(aux->vcf_header, vrec, "DUPLEX_DEPTH", (void *)&n_duplex, 1);
                bcf_update_info(aux->vcf_header, vrec, "BMF_VET", (void *)(&pass_values[0]), vrec->n_allele, BCF_HT_INT);
                bcf_update_info(aux->vcf_header, vrec, "BMF_FAIL", (void *)(&fail_values[0]), vrec->n_allele, BCF_HT_INT);
                bcf_update_info(aux->vcf_header, vrec, "BMF_DUPLEX", (void *)&duplex_values[0], vrec->n_allele, BCF_HT_INT);
                bcf_update_info(aux->vcf_header, vrec, "BMF_UNIOBS", (void *)&uniobs_values[0], vrec->n_allele, BCF_HT_INT);
                bcf_update_info(aux->vcf_header, vrec, "BMF_QUANT", (void *)&quant_est[0], vrec->n_allele, BCF_HT_INT);

                // Pass or fail them individually.
                bcf_write(aux->vcf_ofp, aux->vcf_header, vrec);
                hts_itr_destroy(aux->iter);
            }
            if(vcf_iter) tbx_itr_destroy(vcf_iter);
        }
        LOG_DEBUG("Number skipped for having no reads at position: %i.\n", n_skipped);
        if(bcf_idx) hts_idx_destroy(bcf_idx);
        if(vcf_idx) tbx_destroy(vcf_idx);
        hts_idx_destroy(idx);
        bcf_destroy(vrec);
        return EXIT_SUCCESS;
    }

    int vet_core(vetter_aux_t *aux) {
        return aux->bed ? vet_core_bed(aux): vet_core_nobed(aux);
    }

    int vet_main(int argc, char *argv[])
    {
        if(argc < 3) vetter_usage(EXIT_FAILURE);
        const struct option lopts[] = {
                {"min-family-agreed",         required_argument, nullptr, 'a'},
                {"min-family-size",          required_argument, nullptr, 's'},
                {"min-fraction-agreed",         required_argument, nullptr, 'f'},
                {"min-mapping-quality",         required_argument, nullptr, 'm'},
                {"min-phred-quality",         required_argument, nullptr, 'v'},
                {"min-count",         required_argument, nullptr, 'c'},
                {"min-duplex",         required_argument, nullptr, 'D'},
                {"min-overlap",         required_argument, nullptr, 'O'},
                {"out-vcf",         required_argument, nullptr, 'o'},
                {"bedpath",         required_argument, nullptr, 'b'},
                {"ref",         required_argument, nullptr, 'r'},
                {"padding",         required_argument, nullptr, 'p'},
                {"skip-secondary", no_argument, nullptr, '2'},
                {"skip-supplementary", no_argument, nullptr, 'S'},
                {"skip-qcfail", no_argument, nullptr, 'q'},
                {"skip-improper", no_argument, nullptr, 'P'},
                {"skip-recommended", no_argument, nullptr, 'F'},
                {"max-depth", required_argument, nullptr, 'd'},
                {"emit-bcf", no_argument, nullptr, 'B'},
                {0, 0, 0, 0}
        };
        char vcf_wmode[4] = "w";
        char *outvcf = nullptr, *bed = nullptr;
        int c;
        int padding = 0, output_bcf = 0;
        // Defaults to outputting textual (vcf)
        htsFormat open_fmt = {sequence_data, bam, {1, 3}, gzip, 0, nullptr};
        vetter_aux_t aux = {0};
        aux.min_count = 1;
        aux.max_depth = (1 << 18); // Default max depth

        while ((c = getopt_long(argc, argv, "D:q:r:2:S:d:a:s:m:p:f:b:v:o:O:c:A:BP?hVwF", lopts, nullptr)) >= 0) {
            switch (c) {
            case 'B': output_bcf = 1; break;
            case 'a': aux.minFA = atoi(optarg); break;
            case 'A': NUM_PREALLOCATED_ALLELES = strtoull(optarg, 0, 0); LOG_DEBUG("Num preallocated: %lu\n", NUM_PREALLOCATED_ALLELES); break;
            case 'c': aux.min_count = atoi(optarg); break;
            case 'D': aux.min_duplex = atoi(optarg); break;
            case 's': aux.minFM = atoi(optarg); break;
            case 'm': aux.minmq = atoi(optarg); break;
            case 'v': aux.minPV = atoi(optarg); break;
            case '2': aux.skip_flag |= BAM_FSECONDARY; break;
            case 'S': aux.skip_flag |= BAM_FSUPPLEMENTARY; break;
            case 'q': aux.skip_flag |= BAM_FQCFAIL;break;
            case 'r': aux.skip_flag |= BAM_FDUP; break;
            case 'P': aux.skip_improper = 1; break;
            case 'p': padding = atoi(optarg); break;
            case 'd': aux.max_depth = atoi(optarg); break;
            case 'f': aux.min_fr = (float)atof(optarg); break;
            case 'F': aux.skip_flag |= (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP); break;
            case 'b': bed = strdup(optarg); break;
            case 'o': outvcf = strdup(optarg); break;
            case 'O': aux.min_overlap = atoi(optarg); break;
            case 'V': aux.vet_all = 1; break;
            case 'h': case '?': vetter_usage(EXIT_SUCCESS);
            }
        }

        if(optind + 1 >= argc) vetter_error("Insufficient arguments. Input bam required!\n", EXIT_FAILURE);
        // Check for required tags.
        if(aux.minAF) dlib::check_bam_tag_exit(argv[optind + 1], "AF");
        for(auto tag : {"FA", "FM", "FP", "PV"})
            dlib::check_bam_tag_exit(argv[optind + 1], tag);


        strcpy(vcf_wmode, output_bcf ? "wb": "w");
        if(!outvcf) outvcf = strdup("-");
        if(strcmp(outvcf, "-") == 0) LOG_DEBUG("Emitting to stdout in %s format.\n", output_bcf ? "bcf": "vcf");
        // Open bam
        aux.fp = sam_open_format(argv[optind + 1], "r", &open_fmt);
        if(!aux.fp) LOG_EXIT("Could not open input bam %s. Abort!\n", argv[optind + 1]);
        aux.header = sam_hdr_read(aux.fp);

        // Open input vcf
        if(!aux.header || aux.header->n_targets == 0)
            LOG_EXIT("Could not read header from bam %s. Abort!\n", argv[optind + 1]);
        // Open bed file
        // if no bed provided, do whole genome.
        if(bed) aux.bed = dlib::parse_bed_hash(bed, aux.header, padding);
        //else LOG_EXIT("No bed file provided. Required. Abort!\n");

        if((aux.vcf_fp = vcf_open(argv[optind], "r")) == nullptr) LOG_EXIT("Could not open input vcf (%s).\n", argv[optind]);
        if((aux.vcf_header = bcf_hdr_read(aux.vcf_fp)) == nullptr) LOG_EXIT("Could not read variant header from file (%s).\n", aux.vcf_fp->fn);

        // Add lines to header
        for(unsigned i = 0; i < COUNT_OF(bmf_header_lines); ++i) {
            LOG_DEBUG("Adding header line %s.\n", bmf_header_lines[i]);
            if(bcf_hdr_append(aux.vcf_header, bmf_header_lines[i]))
                LOG_EXIT("Could not add header line '%s'. Abort!\n", bmf_header_lines[i]);
        }
        bcf_hdr_printf(aux.vcf_header, "##bed_filename=\"%s\"", bed ? bed: "FullGenomeAnalysis");
        kstring_t tmpstr = {0};
        ksprintf(&tmpstr, "##cmdline=");
        kputs("bmftools ", &tmpstr);
        for(int i = 0; i < argc; ++i) {
            kputs(argv[i], &tmpstr);
            kputc(' ', &tmpstr);
        }
        bcf_hdr_append(aux.vcf_header, tmpstr.s);
        tmpstr.l = 0;
        // Add in settings
        free(tmpstr.s);
        bcf_hdr_printf(aux.vcf_header, "##bmftools_version=\"%s\"", BMF_VERSION);
        std::string timestring("", 16uL);
        dlib::string_fmt_time(timestring);
        bcf_hdr_printf(aux.vcf_header, "##StartTime=\"%s\"", timestring.c_str());
        dlib::bcf_add_bam_contigs(aux.vcf_header, aux.header);

        // Open output vcf

        if((aux.vcf_ofp = vcf_open(outvcf, vcf_wmode)) == nullptr)
            LOG_EXIT("Could not open output vcf '%s' for writing. Abort!\n", outvcf);
        bcf_hdr_write(aux.vcf_ofp, aux.vcf_header);

        // Open out vcf
        int ret = vet_core(&aux);
        if(ret) fprintf(stderr, "[E:%s:%d] vet_core returned non-zero exit status '%i'. Abort!\n",
                        __func__, __LINE__, ret);
        else LOG_INFO("Successfully completed bmftools vet!\n");
        sam_close(aux.fp);
        bam_hdr_destroy(aux.header);
        vcf_close(aux.vcf_fp);
        vcf_close(aux.vcf_ofp);
        bcf_hdr_destroy(aux.vcf_header);
        if(aux.bed) dlib::bed_destroy_hash(aux.bed);
        cond_free(outvcf);
        cond_free(bed);
        return ret;
    }

} /* namespace bmf */
