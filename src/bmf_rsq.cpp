#include "bmf_rsq.h"
#include <string.h>
#include <getopt.h>
#include "dlib/cstr_util.h"
#include "include/igamc_cephes.h" /// for igamc
#include <algorithm>

namespace BMF {

    struct rsq_aux_t {
        FILE *fqh;
        samFile *in;
        samFile *out;
        int mmlim; // Mismatch failure threshold.
        uint32_t cmpkey:1; // 0 for pos, 1 for unclipped start position
        uint32_t is_se:1;
        uint32_t write_supp:1; // Write reads with supplementary alignments
        uint32_t infer:1; // Use inference instead of barcodes. Fails on already-marked barcoded datasets.
        bam_hdr_t *hdr; // BAM header
        std::unordered_map<std::string, std::string> realign_pairs;
    };

    static const std::function<int (bam1_t *, bam1_t *)> fns[4] = {&same_stack_pos, &same_stack_pos_se,
                                                                   &same_stack_ucs, &same_stack_ucs_se};

    inline void bam2ffq(bam1_t *b, FILE *fp)
    {
        int i;
        uint8_t *rvdata;
        kstring_t ks{0, 120uL, (char *)malloc(120uL)};
        ks.l = 1, ks.s[0] ='@', ks.s[1] = '\0';
        kputs(bam_get_qname(b), &ks);
        kputsnl(" PV:B:I", &ks);
        auto fa((uint32_t *)dlib::array_tag(b, "FA"));
        auto pv((uint32_t *)dlib::array_tag(b, "PV"));
        for(i = 0; i < b->core.l_qseq; ++i) ksprintf(&ks, ",%u", pv[i]);
        kputsnl("\tFA:B:I", &ks);
        for(i = 0; i < b->core.l_qseq; ++i) ksprintf(&ks, ",%u", fa[i]);
        ksprintf(&ks, "\tFM:i:%i\tFP:i:%i", bam_itag(b, "FM"), bam_itag(b, "FP"));
        write_tag_if_found(rvdata, b, "RV", ks);
        write_tag_if_found(rvdata, b, "NC", ks);
        write_tag_if_found(rvdata, b, "DR", ks);
        write_tag_if_found(rvdata, b, "NP", ks);
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
        int n_changed(dlib::int_tag_zero(pdata) + dlib::int_tag_zero(bdata));
        if(pdata) bam_aux_del(p, pdata);
        // If the collapsed observation is now duplex but wasn't before, this updates the DR tag.
        if(pTMP != pFM && pTMP != 0 && (pdata = bam_aux_get(p, "DR")) != nullptr && bam_aux2i(pdata) == 0) {
            pTMP = 1;
            bam_aux_del(p, pdata);
            bam_aux_append(p, "DR", 'i', sizeof(int), (uint8_t *)&pTMP);
        }
        if((pdata = bam_aux_get(p, "NP")) != nullptr) {
            bdata = bam_aux_get(b, "NP");
            pTMP = bam_aux2i(pdata) + (bdata ? bam_aux2i(bdata) : 1);
            bam_aux_del(p, pdata);
        } else {
            bdata = bam_aux_get(b, "NP");
            pTMP = (bdata ? bam_aux2i(bdata) : 1) + 1;
        }
        bam_aux_append(p, "NP", 'i', sizeof(int), reinterpret_cast<uint8_t *>(&pTMP));
        uint32_t *bPV((uint32_t *)dlib::array_tag(b, "PV")); // Length of this should be b->l_qseq
        uint32_t *pPV((uint32_t *)dlib::array_tag(p, "PV"));
        uint32_t *bFA((uint32_t *)dlib::array_tag(b, "FA"));
        uint32_t *pFA((uint32_t *)dlib::array_tag(p, "FA"));
        uint8_t *bSeq(bam_get_seq(b));
        uint8_t *pSeq(bam_get_seq(p));
        uint8_t *bQual(bam_get_qual(b));
        uint8_t *pQual(bam_get_qual(p));
        const int qlen = p->core.l_qseq;
        int8_t ps, bs;

        if(p->core.flag & (BAM_FREVERSE)) {
            int qleni1;
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
                if(pPV[i] > 2) {
                    if((uint32_t)(pQual[i]) > pPV[i]) pQual[i] = (uint8_t)pPV[i];
                } else {
                    pFA[i] = 0;
                    pPV[i] = 0;
                    pQual[i] = 2;
                    n_base(pSeq, i);
                }
            }
        }
        bam_aux_append(p, "NC", 'i', sizeof(int), (uint8_t *)&n_changed);
    }


    void write_stack_se(dlib::tmp_stack_t *stack, rsq_aux_t *settings)
    {
        //size_t n = 0;
        uint8_t *data;
        for(unsigned i = 0; i < stack->n; ++i) {
            if(stack->a[i]) {
                if((data = bam_aux_get(stack->a[i], "NC")) != nullptr) {
                    LOG_DEBUG("Trying to write.\n");
                    bam2ffq(stack->a[i], settings->fqh);
                } else {
                    sam_write1(settings->out, settings->hdr, stack->a[i]);
                }
                bam_destroy1(stack->a[i]), stack->a[i] = nullptr;
            }
        }
    }


    void write_stack(dlib::tmp_stack_t *stack, rsq_aux_t *settings)
    {
        if(settings->is_se) return write_stack_se(stack, settings);
        //size_t n = 0;
        LOG_DEBUG("Starting to write stack\n");
        uint8_t *data;
        for(unsigned i = 0; i < stack->n; ++i) {
            if(stack->a[i]) {
                LOG_DEBUG("Starting to work on this read.\n");
                data = bam_aux_get(stack->a[i], "NC");
                LOG_DEBUG("Got data.\n");
                if(data) {
                    LOG_DEBUG("Trying to write.\n");
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
                    LOG_DEBUG("Trying to write write supp or stuff.\n");
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
                    LOG_DEBUG("About to output bam\n");
                    uint8_t *data;
                    if((data = bam_aux_get(stack->a[i], "MU")) != nullptr)
                        bam_aux_del(stack->a[i], data);
                    sam_write1(settings->out, settings->hdr, stack->a[i]);
                }
                bam_destroy1(stack->a[i]), stack->a[i] = nullptr;
            } else {LOG_DEBUG("DNE\n");}
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


    static inline void flatten_stack_linear(dlib::tmp_stack_t *stack, int mmlim)
    {
        std::sort(stack->a, stack->a + stack->n, [](bam1_t *a, bam1_t *b) {
                return a ? (b ? 0: 1): b ? strcmp(bam_get_qname(a), bam_get_qname(b)): 0;
        });
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
                if(stack->a[i]->core.l_qname != stack->a[j]->core.l_qname)
                    continue;
                if(stringhd(bam_get_qname(stack->a[i]), bam_get_qname(stack->a[j])) < mmlim) {
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

    static inline void flatten_stack_infer(dlib::tmp_stack_t *stack, int mmlim)
    {
        std::sort(stack->a, stack->a + stack->n, [](bam1_t *a, bam1_t *b) {
                return a ? (b ? 0: 1): b ? strcmp(bam_get_qname(a), bam_get_qname(b)): 0;
        });
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

    static const char *sorted_order_strings[2] = {"positional_rescue", "unclipped_rescue"};

    void rsq_core(rsq_aux_t *settings, dlib::tmp_stack_t *stack)
    {
        // This selects the proper function to use for deciding if reads belong in the same stack.
        // It chooses the single-end or paired-end based on is_se and the bmf or pos based on cmpkey.
        const std::function<int (bam1_t *, bam1_t *)> fn = fns[settings->is_se | (settings->cmpkey<<1)];
        if(strcmp(dlib::get_SO(settings->hdr).c_str(), sorted_order_strings[settings->cmpkey]))
            LOG_EXIT("Sort order (%s) is not expected %s for rescue mode. Abort!\n",
                     dlib::get_SO(settings->hdr).c_str(), sorted_order_strings[settings->cmpkey]);
        bam1_t *b = bam_init1();
        uint64_t count = 0;
        while (LIKELY(sam_read1(settings->in, settings->hdr, b) >= 0)) {
            if(UNLIKELY(++count % 1000000 == 0)) LOG_INFO("Records read: %lu.\n", count);
            if(b->core.flag & (BAM_FSUPPLEMENTARY | BAM_FSECONDARY | BAM_FUNMAP | BAM_FMUNMAP)) {
                sam_write1(settings->out, settings->hdr, b);
                continue;
            }
            //LOG_DEBUG("Read a read!\n");
            if(stack->n == 0 || fn(b, *stack->a) == 0) {
                //LOG_DEBUG("Flattening stack\n");
                // New stack -- flatten what we have and write it out.
                flatten_stack_linear(stack, settings->mmlim); // Change this later if the chemistry necessitates it.
                write_stack(stack, settings);
                stack->n = 1;
                stack->a[0] = bam_dup1(b);
                //LOG_DEBUG("Flattened stack\n");
            } else {
                //LOG_DEBUG("Stack now has size %lu.\n", stack->n);
                stack_insert(stack, b);
            }
        }
        flatten_stack_linear(stack, settings->mmlim); // Change this later if the chemistry necessitates it.
        write_stack(stack, settings);
        stack->n = 0;
        bam_destroy1(b);
        // Handle any unpaired reads, though there shouldn't be any in real datasets.
        LOG_DEBUG("Number of orphan reads: %lu.\n", settings->realign_pairs.size());
        if(settings->realign_pairs.size()) {
            LOG_WARNING("There shouldn't be orphan reads in real datasets. Number found: %lu\n", settings->realign_pairs.size());
#if !NDEBUG
            for(auto pair: settings->realign_pairs)
                puts(pair.second.c_str());
#endif
        }
    }

    static const std::vector<uint32_t> ONES(300uL, 1);

    inline void add_dummy_tags(bam1_t *b)
    {
        const int one(1);
        uint8_t *d;
        static std::vector<uint32_t> pvbuf;
        pvbuf.reserve(b->core.l_qseq);
        // Set the read to be a singleton
        if((d = bam_aux_get(b, "FM")) == nullptr)
            bam_aux_append(b, "FM", 'i', sizeof(int), const_cast<uint8_t *>(reinterpret_cast<const uint8_t *>(&one)));
        // Pass the read
        if((d = bam_aux_get(b, "FP")) == nullptr)
            bam_aux_append(b, "FP", 'i', sizeof(int), const_cast<uint8_t *>(reinterpret_cast<const uint8_t *>(&one)));
        if((d = bam_aux_get(b, "FA")) == nullptr) {
            const uint8_t *qual = bam_get_qual(b);
            if(b->core.flag & BAM_FREVERSE) {
                const uint8_t *end = bam_get_qual(b) + b->core.l_qseq;
                while(end >= qual)
                    pvbuf.push_back(static_cast<uint32_t>(*--end));

            } else {
                for(int i = 0; i < b->core.l_qseq; ++i)
                    pvbuf.push_back(static_cast<uint32_t>(qual[i]));
            }
            bam_aux_append(b, "FA", 'B', b->core.l_qseq * sizeof(uint32_t), const_cast<uint8_t *>(reinterpret_cast<const uint8_t *>(ONES.data())));
            bam_aux_append(b, "PV", 'B', b->core.l_qseq * sizeof(uint32_t), const_cast<uint8_t *>(reinterpret_cast<const uint8_t *>(pvbuf.data())));
        }
    }

    void infer_core(rsq_aux_t *settings, dlib::tmp_stack_t *stack)
    {
        // This selects the proper function to use for deciding if reads belong in the same stack.
        // It chooses the single-end or paired-end based on is_se and the bmf or pos based on cmpkey.
        const std::function<int (bam1_t *, bam1_t *)> fn = fns[settings->is_se | (settings->cmpkey<<1)];
        if(strcmp(dlib::get_SO(settings->hdr).c_str(), sorted_order_strings[settings->cmpkey]))
            LOG_EXIT("Sort order (%s) is not expected %s for rescue mode. Abort!\n",
                     dlib::get_SO(settings->hdr).c_str(), sorted_order_strings[settings->cmpkey]);
        bam1_t *b = bam_init1();
        uint64_t count = 0;
        while (LIKELY(sam_read1(settings->in, settings->hdr, b) >= 0)) {
            if(UNLIKELY(++count % 1000000 == 0)) LOG_INFO("Records read: %lu.\n", count);
            add_dummy_tags(b);
            if(b->core.flag & (BAM_FSUPPLEMENTARY | BAM_FSECONDARY | BAM_FUNMAP | BAM_FMUNMAP)) {
                sam_write1(settings->out, settings->hdr, b);
                continue;
            }
            //LOG_DEBUG("Read a read!\n");
            if(stack->n == 0 || fn(b, *stack->a) == 0) {
                LOG_DEBUG("Flattening stack\n");
                // New stack -- flatten what we have and write it out.
                flatten_stack_infer(stack, settings->mmlim); // Change this later if the chemistry necessitates it.
                LOG_DEBUG("Flattened!\n")
                write_stack(stack, settings);
                LOG_DEBUG("Written\n");
                stack->n = 1;
                stack->a[0] = bam_dup1(b);
                //LOG_DEBUG("Flattened stack\n");
            } else {
                LOG_DEBUG("Stack now has size %lu.\n", stack->n);
                stack_insert(stack, b);
            }
        }
        LOG_DEBUG("Loop exit.\n");
        flatten_stack_infer(stack, settings->mmlim); // Change this later if the chemistry necessitates it.
        write_stack(stack, settings);
        stack->n = 0;
        bam_destroy1(b);
        // Handle any unpaired reads, though there shouldn't be any in real datasets.
        LOG_DEBUG("Number of orphan reads: %lu.\n", settings->realign_pairs.size());
        if(settings->realign_pairs.size()) {
            LOG_WARNING("There shouldn't be orphan reads in real datasets. Number found: %lu\n", settings->realign_pairs.size());
#if !NDEBUG
            for(auto pair: settings->realign_pairs)
                puts(pair.second.c_str());
#endif
        }
    }

    void bam_rsq_bookends(rsq_aux_t *settings)
    {
        dlib::tmp_stack_t stack{0, STACK_START, (bam1_t **)malloc(STACK_START * sizeof(bam1_t *))};
        if(settings->infer) infer_core(settings, &stack);
        else rsq_core(settings, &stack);
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
                        "-i      Flag to ignore barcodes and infer solely by positional information. [THIS IS BROKEN]\n"
                        "This flag adds artificial auxiliary tags to treat unbarcoded reads as if they were singletons.\n"
                );
        return retcode;
    }


    int rsq_main(int argc, char *argv[])
    {
        int c;
        char wmode[4] = "wb";

        rsq_aux_t settings = {0};
        settings.mmlim = 2;
        settings.cmpkey = cmpkey::POSITION;

        char *fqname = nullptr;

        if(argc < 3) return rsq_usage(EXIT_FAILURE);

        while ((c = getopt(argc, argv, "l:f:t:iSHush?")) >= 0) {
            switch (c) {
            case 's': settings.write_supp = 1; break;
            case 'S': settings.is_se = 1; break;
            case 'u':
                settings.cmpkey = cmpkey::UNCLIPPED;
                LOG_INFO("Unclipped start position chosen for cmpkey.\n");
                break;
            case 't': settings.mmlim = atoi(optarg); break;
            case 'f': fqname = optarg; break;
            case 'l': wmode[2] = atoi(optarg)%10 + '0';break;
            case 'i': settings.infer = 1; break;
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


        if(settings.cmpkey == cmpkey::UNCLIPPED)
            if(!settings.is_se)
                dlib::check_bam_tag_exit(argv[optind], "MU");
        if(!settings.infer)
            for(const char *tag: {"FM", "FA", "PV", "FP"})
                dlib::check_bam_tag_exit(argv[optind], tag);
        settings.in = sam_open(argv[optind], "r");
        settings.hdr = sam_hdr_read(settings.in);

        if (settings.hdr == nullptr || settings.hdr->n_targets == 0)
            LOG_EXIT("input SAM does not have header. Abort!\n");

        dlib::add_pg_line(settings.hdr, argc, argv, "bmftools rsq", BMF_VERSION, "bmftools",
                "Uses positional information to rescue reads with errors in the barcode.");
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
