#include "bmf_rsq.h"

namespace BMF {

    struct rsq_aux_t {
        FILE *fqh;
        samFile *in;
        samFile *out;
        int mmlim; // Mismatch failure threshold.
        uint32_t cmpkey:1; // 0 for pos, 1 for unclipped start position
        uint32_t is_se:1;
        uint32_t write_supp:1; // Write reads with supplementary alignments
        bam_hdr_t *hdr; // BAM header
        std::unordered_map<std::string, std::string> realign_pairs;
    };

    static const std::function<int (bam1_t *, bam1_t *)> fns[4] = {&same_stack_pos, &same_stack_pos_se,
                                                                   &same_stack_ucs, &same_stack_ucs_se};

    inline void bam2ffq(bam1_t *b, FILE *fp)
    {
        int i;
        uint8_t *rvdata;
        kstring_t ks{0, 0, nullptr};
        kputc('@', &ks);
        kputs(bam_get_qname(b), &ks);
        kputs(" PV:B:I", &ks);
        auto fa((uint32_t *)dlib::array_tag(b, "FA"));
        auto pv((uint32_t *)dlib::array_tag(b, "PV"));
        for(i = 0; i < b->core.l_qseq; ++i) ksprintf(&ks, ",%u", pv[i]);
        kputs("\tFA:B:I", &ks);
        for(i = 0; i < b->core.l_qseq; ++i) ksprintf(&ks, ",%u", fa[i]);
        ksprintf(&ks, "\tFM:i:%i\tFP:i:%i",
                 bam_itag(b, "FM"), bam_itag(b, "FP"));
        if((rvdata = bam_aux_get(b, "RV")) != nullptr)
            ksprintf(&ks, "\tRV:i:%i", bam_aux2i(rvdata));
        if((rvdata = bam_aux_get(b, "NC")) != nullptr)
            ksprintf(&ks, "\tNC:i:%i", bam_aux2i(rvdata));
        if((rvdata = bam_aux_get(b, "DR")) != nullptr)
            ksprintf(&ks, "\tDR:i:%i", bam_aux2i(rvdata));
        kputc('\n', &ks);
        uint8_t *seq(bam_get_seq(b));
        char *seqbuf((char *)malloc(b->core.l_qseq + 1));
        for (i = 0; i < b->core.l_qseq; ++i) seqbuf[i] = seq_nt16_str[bam_seqi(seq, i)];
        seqbuf[i] = '\0';
        if (b->core.flag & BAM_FREVERSE) { // reverse complement
            for(i = 0; i < b->core.l_qseq>>1; ++i) {
                const int8_t t = seqbuf[b->core.l_qseq - i - 1];
                seqbuf[b->core.l_qseq - i - 1] = nuc_cmpl(seqbuf[i]);
                seqbuf[i] = nuc_cmpl(t);
            }
            if(b->core.l_qseq&1) seqbuf[i] = nuc_cmpl(seqbuf[i]);
        }
        seqbuf[b->core.l_qseq] = '\0';
        assert(strlen(seqbuf) == (uint64_t)b->core.l_qseq);
        kputs(seqbuf, &ks);
        kputs("\n+\n", &ks);
        uint8_t *qual(bam_get_qual(b));
        for(i = 0; i < b->core.l_qseq; ++i) seqbuf[i] = 33 + qual[i];
        if (b->core.flag & BAM_FREVERSE) { // reverse
            for (i = 0; i < b->core.l_qseq>>1; ++i) {
                const int8_t t = seqbuf[b->core.l_qseq - 1 - i];
                seqbuf[b->core.l_qseq - 1 - i] = seqbuf[i];
                seqbuf[i] = t;
            }
        }
        kputs(seqbuf, &ks), free(seqbuf);
        kputc('\n', &ks);
        fputs(ks.s, fp), free(ks.s);
    }

    inline int switch_names(char *n1, char *n2) {
        for(;*n1;++n1, ++n2)
            if(*n1 != *n2)
                return *n1 < *n2;
        return 0; // If identical, don't switch. Should never happen.
    }

    inline void bam2ffq(bam1_t *b, FILE *fp, std::string& qname)
    {
        int i;
        uint8_t *rvdata;
        kstring_t ks{0, 0, nullptr};
        kputc('@', &ks);
        kputs(qname.c_str(), &ks);
        kputs(" PV:B:I", &ks);
        auto fa((uint32_t *)dlib::array_tag(b, "FA"));
        auto pv((uint32_t *)dlib::array_tag(b, "PV"));
        for(i = 0; i < b->core.l_qseq; ++i) ksprintf(&ks, ",%u", pv[i]);
        kputs("\tFA:B:I", &ks);
        for(i = 0; i < b->core.l_qseq; ++i) ksprintf(&ks, ",%u", fa[i]);
        ksprintf(&ks, "\tFM:i:%i\tFP:i:%i\tNC:i:%i",
                bam_itag(b, "FM"), bam_itag(b, "FP"), bam_itag(b, "NC"));
        if((rvdata = bam_aux_get(b, "RV")) != nullptr)
            ksprintf(&ks, "\tRV:i:%i", bam_aux2i(rvdata));
        kputc('\n', &ks);
        uint8_t *seq(bam_get_seq(b));
        char *seqbuf((char *)malloc(b->core.l_qseq + 1));
        for (i = 0; i < b->core.l_qseq; ++i) seqbuf[i] = seq_nt16_str[bam_seqi(seq, i)];
        seqbuf[i] = '\0';
        if (b->core.flag & BAM_FREVERSE) { // reverse complement
            for(i = 0; i < b->core.l_qseq>>1; ++i) {
                const int8_t t = seqbuf[b->core.l_qseq - i - 1];
                seqbuf[b->core.l_qseq - i - 1] = nuc_cmpl(seqbuf[i]);
                seqbuf[i] = nuc_cmpl(t);
            }
            if(b->core.l_qseq&1) seqbuf[i] = nuc_cmpl(seqbuf[i]);
        }
        seqbuf[b->core.l_qseq] = '\0';
        assert(strlen(seqbuf) == (uint64_t)b->core.l_qseq);
        kputs(seqbuf, &ks);
        kputs("\n+\n", &ks);
        uint8_t *qual(bam_get_qual(b));
        for(i = 0; i < b->core.l_qseq; ++i) seqbuf[i] = 33 + qual[i];
        if (b->core.flag & BAM_FREVERSE) { // reverse
            for (i = 0; i < b->core.l_qseq>>1; ++i) {
                const int8_t t = seqbuf[b->core.l_qseq - 1 - i];
                seqbuf[b->core.l_qseq - 1 - i] = seqbuf[i];
                seqbuf[i] = t;
            }
        }
        kputs(seqbuf, &ks), free(seqbuf);
        kputc('\n', &ks);
        fputs(ks.s, fp), free(ks.s);
    }

    void update_bam1(bam1_t *p, bam1_t *b)
    {
        uint8_t *bdata(bam_aux_get(b, "FM"));
        uint8_t *pdata(bam_aux_get(p, "FM"));
        if(UNLIKELY(!bdata || !pdata)) {
            fprintf(stderr, "Required FM tag not found. Abort mission!\n");
            exit(EXIT_FAILURE);
        }
        int bFM(bam_aux2i(bdata));
        int pFM(bam_aux2i(pdata));
        int pTMP;
        if(switch_names(bam_get_qname(p), bam_get_qname(b))) {
            memcpy(bam_get_qname(p), bam_get_qname(b), b->core.l_qname);
            assert(strlen(bam_get_qname(p)) == strlen(bam_get_qname(b)));
        }
        pFM += bFM;
        bam_aux_del(p, pdata);
        bam_aux_append(p, "FM", 'i', sizeof(int), (uint8_t *)&pFM);
        if((pdata = bam_aux_get(p, "RV")) != nullptr) {
            pTMP = bam_aux2i(pdata) + bam_itag(b, "RV");
            bam_aux_del(p, pdata);
            bam_aux_append(p, "RV", 'i', sizeof(int), (uint8_t *)&pTMP);
        }
        // Handle NC (Number Changed) tag
        pdata = bam_aux_get(p, "NC");
        bdata = bam_aux_get(b, "NC");
        int n_changed{dlib::int_tag_zero(pdata) + dlib::int_tag_zero(bdata)};
        if(pdata) bam_aux_del(p, pdata);
        pdata = bam_aux_get(p, "DR");
        if(bam_aux2i(pdata) == 0 && pTMP != pFM && pTMP != 0) {
            pTMP = 1;
            bam_aux_del(p, pdata);
            bam_aux_append(p, "DR", 'i', sizeof(int), (uint8_t *)&pTMP);
        }
        uint32_t *bPV((uint32_t *)dlib::array_tag(b, "PV")); // Length of this should be b->l_qseq
        uint32_t *pPV((uint32_t *)dlib::array_tag(p, "PV"));
        uint32_t *bFA((uint32_t *)dlib::array_tag(b, "FA")); // Length of this should be b->l_qseq
        uint32_t *pFA((uint32_t *)dlib::array_tag(p, "FA"));
        uint8_t *bSeq(bam_get_seq(b));
        uint8_t *pSeq(bam_get_seq(p));
        uint8_t *bQual(bam_get_qual(b));
        uint8_t *pQual(bam_get_qual(p));
        const int qlen = p->core.l_qseq;

        if(p->core.flag & (BAM_FREVERSE)) {
            int qleni1;
            int8_t ps, bs;
            for(int i = 0; i < qlen; ++i) {
                qleni1 = qlen - i - 1;
                ps = bam_seqi(pSeq, qleni1);
                bs = bam_seqi(bSeq, qleni1);
                if(ps == bs) {
                    pPV[i] = agreed_pvalues(pPV[i], bPV[i]);
                    pFA[i] += bFA[i];
                    if(bQual[qleni1] > pQual[qleni1]) pQual[qleni1] = bQual[qleni1];
                } else if(ps == dlib::htseq::HTS_N) {
                    bam_set_base(pSeq, bSeq, qleni1);
                    pFA[i] = bFA[i];
                    pPV[i] = bPV[i];
                    pQual[qleni1] = bQual[qleni1];
                    ++n_changed; // Note: goes from N to a useable nucleotide.
                    continue;
                } else if(bs == dlib::htseq::HTS_N) {
                    continue;
                } else {
                    if(pPV[i] > bPV[i]) {
                        bam_set_base(pSeq, bSeq, qleni1);
                        pPV[i] = disc_pvalues(pPV[i], bPV[i]);
                    } else pPV[i] = disc_pvalues(bPV[i], pPV[i]);
                    pFA[i] = bFA[i];
                    pQual[qleni1] = bQual[qleni1];
                    ++n_changed;
                }
                if(pPV[i] < 3) {
                    pFA[i] = 0;
                    pPV[i] = 0;
                    pQual[qleni1] = 2;
                    n_base(pSeq, qleni1);
                    continue;
                }
                if((uint32_t)(pQual[qleni1]) > pPV[i]) pQual[qleni1] = (uint8_t)pPV[i];
            }
        } else {
            int8_t ps, bs;
            for(int i = 0; i < qlen; ++i) {
                ps = bam_seqi(pSeq, i);
                bs = bam_seqi(bSeq, i);
                if(ps == bs) {
                    pPV[i] = agreed_pvalues(pPV[i], bPV[i]);
                    pFA[i] += bFA[i];
                    if(bQual[i] > pQual[i]) pQual[i] = bQual[i];
                } else if(ps == dlib::htseq::HTS_N) {
                    bam_set_base(pSeq, bSeq, i);
                    pFA[i] = bFA[i];
                    pPV[i] = bPV[i];
                    ++n_changed; // Note: goes from N to a useable nucleotide.
                    continue;
                } else if(bs == dlib::htseq::HTS_N) {
                    continue;
                } else {
                    if(pPV[i] > bPV[i]) {
                        bam_set_base(pSeq, bSeq, i);
                        pPV[i] = disc_pvalues(pPV[i], bPV[i]);
                    } else pPV[i] = disc_pvalues(bPV[i], pPV[i]);
                    pFA[i] = bFA[i];
                    pQual[i] = bQual[i];
                    ++n_changed;
                }
                if(pPV[i] < 3) {
                    pFA[i] = 0;
                    pPV[i] = 0;
                    pQual[i] = 2;
                    n_base(pSeq, i);
                    continue;
                } else if((uint32_t)(pQual[i]) > pPV[i]) pQual[i] = (uint8_t)pPV[i];
            }
        }
        bam_aux_append(p, "NC", 'i', sizeof(int), (uint8_t *)&n_changed);
    }


    void write_stack_se(dlib::tmp_stack_t *stack, rsq_aux_t *settings);

    void write_stack(dlib::tmp_stack_t *stack, rsq_aux_t *settings)
    {
        if(settings->is_se) return write_stack_se(stack, settings);
        //size_t n = 0;
        uint8_t *data;
        for(unsigned i = 0; i < stack->n; ++i) {
            if(stack->a[i]) {
                if((data = bam_aux_get(stack->a[i], "NC")) != nullptr) {
                    std::string&& qname = bam_get_qname(stack->a[i]);
                    if(settings->realign_pairs.find(qname) == settings->realign_pairs.end()) {
                        settings->realign_pairs[qname] = dlib::bam2cppstr(stack->a[i]);
                    } else {
                        // Make sure the read names/barcodes match.
                        // Write read 1 out first.
                        if(stack->a[i]->core.flag & BAM_FREAD2) {
                            fputs(settings->realign_pairs[qname].c_str(), settings->fqh);
                            bam2ffq(stack->a[i], settings->fqh);
                        } else {
                            bam2ffq(stack->a[i], settings->fqh);
                            fputs(settings->realign_pairs[qname].c_str(), settings->fqh);
                        }
                        // Clear entry, as there can only be two.
                        settings->realign_pairs.erase(qname);
                    }
                } else if(settings->write_supp && (bam_aux_get(stack->a[i], "SA") || bam_aux_get(stack->a[i], "ms"))) {
                    // Has an SA or ms tag, meaning that the read or its mate had a supplementary alignment
                    std::string&& qname = bam_get_qname(stack->a[i]);
                    if(settings->realign_pairs.find(qname) == settings->realign_pairs.end()) {
                        settings->realign_pairs[qname] = dlib::bam2cppstr(stack->a[i]);
                    } else {
                        // Make sure the read names/barcodes match.
                        //assert(memcmp(settings->realign_pairs[qname].c_str() + 1, qname.c_str(), qname.size() - 1) == 0);
                        // Write read 1 out first.
                        if(stack->a[i]->core.flag & BAM_FREAD2) {
                            fputs(settings->realign_pairs[qname].c_str(), settings->fqh);
                            //bam2ffq(stack->a[i], settings->fqh, qname);
                            bam2ffq(stack->a[i], settings->fqh);
                        } else {
                            //bam2ffq(stack->a[i], settings->fqh, qname);
                            bam2ffq(stack->a[i], settings->fqh);
                            fputs(settings->realign_pairs[qname].c_str(), settings->fqh);
                        }
                        // Clear entry, as there can only be two.
                        settings->realign_pairs.erase(qname);
                    }
                } else {
                    uint8_t *data;
                    data = bam_aux_get(stack->a[i], "MU");
                    bam_aux_del(stack->a[i], data);
                    sam_write1(settings->out, settings->hdr, stack->a[i]);
                }
                bam_destroy1(stack->a[i]), stack->a[i] = nullptr;
            }
        }
    }


    void write_stack_se(dlib::tmp_stack_t *stack, rsq_aux_t *settings)
    {
        //size_t n = 0;
        uint8_t *data;
        for(unsigned i = 0; i < stack->n; ++i) {
            if(stack->a[i]) {
                if((data = bam_aux_get(stack->a[i], "NC")) != nullptr) {
                    bam2ffq(stack->a[i], settings->fqh);
                } else {
                    sam_write1(settings->out, settings->hdr, stack->a[i]);
                }
                bam_destroy1(stack->a[i]), stack->a[i] = nullptr;
            }
        }
    }

    inline int stringhd(char *a, char *b) {
        int hd = 0;
        while(*a) hd += (*a++ != *b++);
        return hd;
    }

    inline int string_linear(char *a, char *b, int mmlim) {
        int hd = 0;
        while(*a)
            if(*a++ != *b++)
                if(++hd > mmlim)
                    return 0;
        return 1;
    }
    /*
     * Returns true if hd <= mmlim, 0 otherwise.
     */
#define hd_test(a, b, mmlim) stringhd(bam_get_qname(a), bam_get_qname(b)) < mmlim)

#if !NDEBUG
    namespace {
        KHASH_MAP_INIT_INT(hd, uint64_t)
    }

#endif

#if !NDEBUG
    static inline void flatten_stack_linear(dlib::tmp_stack_t *stack, int mmlim, khash_t(hd) *readhds)
#else
    static inline void flatten_stack_linear(dlib::tmp_stack_t *stack, int mmlim)
#endif
    {
        std::sort(stack->a, stack->a + stack->n, [](bam1_t *a, bam1_t *b) {
                return a ? (b ? 0: 1): b ? strcmp(bam_get_qname(a), bam_get_qname(b)): 0;
        });
#if 0
        const uint64_t key = stack->a[0] ? ucs_sort_core_key(stack->a[0]) : -1;
        const int tid = stack->a[0]->core.tid;
        LOG_DEBUG("Get cs\n");
        const int ucs = dlib::get_unclipped_start(stack->a[0]);
        LOG_DEBUG("Get mucs\n");
        const int mucs = bam_itag(stack->a[0], "MU");
        //LOG_DEBUG("%lu", key);
#endif
        for(unsigned i = 0; i < stack->n; ++i) {
            for(unsigned j = i + 1; j < stack->n; ++j) {
                //assert(key == ucs_sort_core_key(stack->a[j]));
                //assert(ucs == dlib::get_unclipped_start(stack->a[j]));
                //assert(tid == stack->a[j]->core.tid);
                //assert(mucs == bam_itag(stack->a[j], "MU"));
                assert(stack->a[i]);
                assert(stack->a[j]);
                if(stack->a[i]->core.l_qseq != stack->a[j]->core.l_qseq)
                    continue;
                if(stringhd(bam_get_qname(stack->a[i]), bam_get_qname(stack->a[j])) < mmlim) {
#if 0
                    const int readhd = read_hd(stack->a[i], stack->a[j]);
                    if(readhd > 10) {
                        char buf1[200];
                        char buf2[200];
                        dlib::bam_seq_cpy(buf1, stack->a[i]);
                        dlib::bam_seq_cpy(buf2, stack->a[j]);
                        LOG_DEBUG("Encountered large read hd %i. Seq1: %s. Seq2: %s. Name1: %s. Name2: %s. Name HD: %i\n",
                                  readhd, buf1, buf2, bam_get_qname(stack->a[i]), bam_get_qname(stack->a[j]), stringhd(bam_get_qname(stack->a[i]), bam_get_qname(stack->a[j])));
                        assert(stringhd(bam_get_qname(stack->a[i]), bam_get_qname(stack->a[j])) < mmlim);
                    }
                    khiter_t k = kh_get(hd, readhds, readhd);
                    int khr;
                    if(k == kh_end(readhds)) {
                        LOG_DEBUG("Encountered new hd for reads: %i.\n", readhd);
                        k = kh_put(hd, readhds, readhd, &khr);
                        kh_val(readhds, k) = 1;
                    } else ++kh_val(readhds, k);
#endif
                    assert(stringhd(bam_get_qname(stack->a[i]), bam_get_qname(stack->a[j])) <= mmlim);
                    //LOG_DEBUG("Flattening %s into %s.\n", bam_get_qname(stack->a[i]), bam_get_qname(stack->a[j]));
                    update_bam1(stack->a[j], stack->a[i]);
                    bam_destroy1(stack->a[i]);
                    stack->a[i] = nullptr;
                    break;
                    // "break" in case there are multiple within hamming distance.
                    // Otherwise, I'll end up having memory mistakes.
                    // Besides, that read set will get merged into the later read in the set.
                }
            }
        }
    }

    static const char *sorted_order_strings[2] = {"positional_rescue", "unclipped_rescue"};

    void rsq_core(rsq_aux_t *settings, dlib::tmp_stack_t *stack)
    {
#if !NDEBUG
        khash_t(hd) *read_hds = kh_init(hd);
#endif
        // This selects the proper function to use for deciding if reads belong in the same stack.
        // It chooses the single-end or paired-end based on is_se and the bmf or pos based on cmpkey.
        const std::function<int (bam1_t *, bam1_t *)> fn = fns[settings->is_se | (settings->cmpkey<<1)];
        if(strcmp(dlib::get_SO(settings->hdr).c_str(), sorted_order_strings[settings->cmpkey]))
            LOG_EXIT("Sort order (%s) is not expected %s for rescue mode. Abort!\n",
                     dlib::get_SO(settings->hdr).c_str(), sorted_order_strings[settings->cmpkey]);
        bam1_t *b = bam_init1();
        // Start stack
        while (LIKELY(sam_read1(settings->in, settings->hdr, b) >= 0)) {
            if(b->core.flag & (BAM_FSUPPLEMENTARY | BAM_FSECONDARY | BAM_FUNMAP | BAM_FMUNMAP)) {
                sam_write1(settings->out, settings->hdr, b);
                continue;
            }
            //LOG_DEBUG("Read a read!\n");
            if(stack->n == 0 || fn(b, *stack->a) == 0) {
                //LOG_DEBUG("Flattening stack\n");
                // New stack -- flatten what we have and write it out.
#if !NDEBUG
                flatten_stack_linear(stack, settings->mmlim, read_hds); // Change this later if the chemistry necessitates it.
#else
                flatten_stack_linear(stack, settings->mmlim); // Change this later if the chemistry necessitates it.
#endif
                write_stack(stack, settings);
                stack->n = 1;
                stack->a[0] = bam_dup1(b);
                //LOG_DEBUG("Flattened stack\n");
            } else {
                //LOG_DEBUG("Stack now has size %lu.\n", stack->n);
                stack_insert(stack, b);
            }
        }
#if !NDEBUG
        flatten_stack_linear(stack, settings->mmlim, read_hds); // Change this later if the chemistry necessitates it.
        for(khiter_t ki = kh_begin(read_hds); ki != kh_end(read_hds); ++ki)
            if(kh_exist(read_hds, ki))
                fprintf(stdout, "%i:%lu\n", kh_key(read_hds, ki), kh_val(read_hds, ki));
        kh_destroy(hd, read_hds);
#else
        flatten_stack_linear(stack, settings->mmlim); // Change this later if the chemistry necessitates it.
#endif
        write_stack(stack, settings);
        stack->n = 0;
        bam_destroy1(b);
        // Handle any unpaired reads, though there shouldn't be any in real datasets.
        LOG_DEBUG("Number of orphan reads: %lu.\n", settings->realign_pairs.size());
        for(auto pair: settings->realign_pairs)
            fputs(pair.second.c_str(), settings->fqh);
    }

    void bam_rsq_bookends(rsq_aux_t *settings)
    {
        LOG_DEBUG("Starting stack!\n");
        dlib::tmp_stack_t stack = {0, STACK_START, (bam1_t **)malloc(STACK_START * sizeof(bam1_t *))};
        rsq_core(settings, &stack);
        free(stack.a);
    }

    int rsq_usage(int retcode)
    {
        fprintf(stderr,
                        "Positional rescue. \n"
                        "Reads with the same start position (or unclipped start, if -u is set) are compared.\n"
                        "If their barcodes are sufficiently similar, they are treated as having originated"
                        "from the same original template molecule.\n"
                        "Usage:  bmftools rsq <input.srt.bam> <output.bam>\n\n"
                        "Flags:\n"
                        "-f      Path for the fastq for reads that need to be realigned. REQUIRED.\n"
                        "-s      Flag to write reads with supplementary alignments . Default: False.\n"
                        "-S      Flag to indicate that this rescue is for single-end data.\n"
                        "-t      Mismatch limit. Default: 2\n"
                        "-l      Set bam compression level. Valid: 0-9. (0 == uncompresed)\n"
                        "-u      Flag to use unclipped start positions instead of pos/mpos for identifying potential duplicates.\n"
                        "Note: -u requires pre-processing with bmftools mark_unclipped.\n"
                );
        return retcode;
    }


    int rsq_main(int argc, char *argv[])
    {
        int c;
        char wmode[4] = "wb";

        rsq_aux_t settings = {0};
        settings.mmlim = 2;
        settings.cmpkey = POSITION;

        char *fqname = nullptr;

        if(argc < 3) return rsq_usage(EXIT_FAILURE);

        while ((c = getopt(argc, argv, "l:f:t:SHush?")) >= 0) {
            switch (c) {
            case 's': settings.write_supp = 1; break;
            case 'S': settings.is_se = 1; break;
            case 'u':
                settings.cmpkey = UNCLIPPED;
                LOG_INFO("Unclipped start position chosen for cmpkey.\n");
                break;
            case 't': settings.mmlim = atoi(optarg); break;
            case 'f': fqname = optarg; break;
            case 'l': wmode[2] = atoi(optarg)%10 + '0';break;
            case '?': case 'h': case 'H': return rsq_usage(EXIT_SUCCESS);
            }
        }
        if (optind + 2 > argc)
            return rsq_usage(EXIT_FAILURE);

        if(!fqname) {
            fprintf(stderr, "Fastq path for rescued reads required. Abort!\n");
            return rsq_usage(EXIT_FAILURE);
        }

        settings.fqh = fopen(fqname, "w");

        if(!settings.fqh)
            LOG_EXIT("Failed to open output fastq for writing. Abort!\n");


        if(settings.cmpkey == UNCLIPPED)
            for(const char *tag: {"MU"})
                dlib::check_bam_tag_exit(argv[optind], tag);
        for(const char *tag: {"FM", "FA", "PV", "FP", "RV"})
            dlib::check_bam_tag_exit(argv[optind], tag);
        settings.in = sam_open(argv[optind], "r");
        settings.hdr = sam_hdr_read(settings.in);

        if (settings.hdr == nullptr || settings.hdr->n_targets == 0)
            LOG_EXIT("input SAM does not have header. Abort!\n");

        LOG_DEBUG("Write mode: %s.\n", wmode);
        settings.out = sam_open(argv[optind+1], wmode);
        if (settings.in == 0 || settings.out == 0)
            LOG_EXIT("fail to read/write input files\n");
        sam_hdr_write(settings.out, settings.hdr);

        bam_rsq_bookends(&settings);
        bam_hdr_destroy(settings.hdr);
        sam_close(settings.in); sam_close(settings.out);
        if(settings.fqh) fclose(settings.fqh);
        LOG_INFO("Successfully completed bmftools rsq.\n");
        return EXIT_SUCCESS;
    }

}
