#include "bmf_rsq.h"

namespace BMF {

    namespace {
        inline std::string make_name(bam1_t *b, size_t n) {
            std::string ret;
            ret.resize(70uL);
            uint8_t *data_su{bam_aux_get(b, "SU")};
            uint8_t *data_mu{bam_aux_get(b, "MU")};
            if(b->core.flag & BAM_FREAD1) {
                stringprintf(ret, "collapsed:%i:%i:%i:%i:%i:%i:%i:%i:%lu",
                             data_su ? bam_aux2i(data_su): b->core.pos,
                             data_mu ? bam_aux2i(data_mu): b->core.mpos, // Unclipped starts for self and mate
                             b->core.tid,
                             b->core.mtid, // Contigs
                             (!!(b->core.flag & (BAM_FREVERSE))),
                             (!!(b->core.flag & (BAM_FMREVERSE))), // Strandedness combinations
                             b->core.l_qseq,
                             bam_itag(b, "LM"),
                             n // Read length of self and mate.
                             );
            } else {
                stringprintf(ret, "collapsed:%i:%i:%i:%i:%i:%i:%i:%i:%lu",
                             data_mu ? bam_aux2i(data_mu): b->core.mpos, // Unclipped starts for self and mate
                             data_su ? bam_aux2i(data_su): b->core.pos,
                             b->core.mtid,
                             b->core.tid,
                             (!!(b->core.flag & (BAM_FMREVERSE))), // Strandedness combinations
                             (!!(b->core.flag & (BAM_FREVERSE))),
                             bam_itag(b, "LM"),
                             b->core.l_qseq,
                             n // Read length of self and mate.
                             );
            }
            return ret;
        }
    }

    struct rsq_aux_t {
        FILE *fqh;
        samFile *in;
        samFile *out;
        int cmpkey; // 0 for pos, 1 for unclipped start position
        int mmlim; // Mismatch failure threshold.
        int is_se;
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
        auto fa{(uint32_t *)dlib::array_tag(b, "FA")};
        auto pv{(uint32_t *)dlib::array_tag(b, "PV")};
        fa = (uint32_t *)dlib::array_tag(b, "FA");
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
        auto fa{(uint32_t *)dlib::array_tag(b, "FA")};
        auto pv{(uint32_t *)dlib::array_tag(b, "PV")};
        fa = (uint32_t *)dlib::array_tag(b, "FA");
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
        if(switch_names(bam_get_qname(p), bam_get_qname(b))) {
            memcpy(bam_get_qname(p), bam_get_qname(b), b->core.l_qname);
            assert(strlen(bam_get_qname(p)) == strlen(bam_get_qname(b)));
        }
        pFM += bFM;
        bam_aux_del(p, pdata);
        bam_aux_append(p, "FM", 'i', sizeof(int), (uint8_t *)&pFM);
        if((pdata = bam_aux_get(p, "RV")) != nullptr) {
            const int pRV = bam_aux2i(pdata) + bam_itag(b, "RV");
            bam_aux_del(p, pdata);
            bam_aux_append(p, "RV", 'i', sizeof(int), (uint8_t *)&pRV);
        }
        // Handle NC (Number Changed) tag
        pdata = bam_aux_get(p, "NC");
        bdata = bam_aux_get(b, "NC");
        int n_changed{dlib::int_tag_zero(pdata) + dlib::int_tag_zero(bdata)};
        if(pdata) bam_aux_del(p, pdata);

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
            for(int i = 0; i < qlen; ++i) {
                const int8_t ps(bam_seqi(pSeq, i));
                const int8_t bs(bam_seqi(bSeq, i));
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


    void write_stack(dlib::tmp_stack_t *stack, rsq_aux_t *settings)
    {
        size_t n = 0;
        for(unsigned i = 0; i < stack->n; ++i) {
            if(stack->a[i]) {
                uint8_t *data;
                if((data = bam_aux_get(stack->a[i], "NC")) != nullptr) {
                    std::string&& qname = make_name(stack->a[i], ++n);
                    LOG_DEBUG("Read name: \t%s\tGenerated name: \t%s\t.\n", bam_get_qname(stack->a[i]), qname.c_str());
                    //LOG_DEBUG("Qname: %s.\n", qname.c_str());
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
                } else sam_write1(settings->out, settings->hdr, stack->a[i]);
                bam_destroy1(stack->a[i]), stack->a[i] = nullptr;
            }
        }
    }

    /*
     * Returns true if hd <= mmlim, 0 otherwise.
     */
    static inline int hd_linear(bam1_t *a, bam1_t *b, int mmlim)
    {
        char *aname = (char *)bam_get_qname(a);
        char *bname = (char *)bam_get_qname(b);
        int hd = 0;
        while(*aname)
            if(*aname++ != *bname++)
                if(++hd > mmlim)
                    return 0;
        return 1;
    }

    static inline void flatten_stack_linear(dlib::tmp_stack_t *stack, rsq_aux_t *settings)
    {
        // Sort by read names to make sure that any progressive rescuing ends at the same name.
        std::sort(stack->a, &stack->a[stack->n], [](const bam1_t *a, const bam1_t *b) {
            if(a) {
                return b ? (int)(strcmp(bam_get_qname(a), bam_get_qname(b)) < 0): 0;
            } else {
                return b ? 1: 0;
            }
            // return a ? (b ? (int)(strcmp(bam_get_qname(a), bam_get_qname(b)) < 0): 0): (b ? 1: 0);
            // Returns 0 if comparing two nulls, and returns true that a nullptr lt a valued name
            // Compares strings otherwise.
        });
        for(unsigned i = 0; i < stack->n; ++i) {
            for(unsigned j = i + 1; j < stack->n; ++j) {
                if(hd_linear(stack->a[i], stack->a[j], settings->mmlim) &&
                        stack->a[i]->core.l_qseq == stack->a[j]->core.l_qseq) {
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
        // This selects the proper function to use for deciding if reads belong in the same stack.
        // It chooses the single-end or paired-end based on is_se and the bmf or pos based on cmpkey.
        std::function<int (bam1_t *, bam1_t *)> fn = fns[settings->is_se | (settings->cmpkey<<1)];
        if(strcmp(dlib::get_SO(settings->hdr).c_str(), sorted_order_strings[settings->cmpkey]))
            LOG_EXIT("Sort order (%s) is not expected %s for rescue mode. Abort!\n",
                     dlib::get_SO(settings->hdr).c_str(), sorted_order_strings[settings->cmpkey]);
        bam1_t *b = bam_init1();
        if(sam_read1(settings->in, settings->hdr, b) < 0)
            LOG_EXIT("Failed to read first record in bam file. Abort!\n");
        // Zoom ahead to first primary alignment in bam.
        while(b->core.flag & (BAM_FSUPPLEMENTARY | BAM_FSECONDARY)) {
            sam_write1(settings->out, settings->hdr, b);
            if(sam_read1(settings->in, settings->hdr, b))
                LOG_EXIT("Could not read first primary alignment in bam (%s). Abort!\n", settings->in->fn);
        }
        // Start stack
        stack_insert(stack, b);
        while (LIKELY(sam_read1(settings->in, settings->hdr, b) >= 0)) {
            if(b->core.flag & (BAM_FSUPPLEMENTARY | BAM_FSECONDARY)) {
                sam_write1(settings->out, settings->hdr, b);
                continue;
            }
            if(fn(b, *stack->a) == 0) {
                // New stack -- flatten what we have and write it out.
                flatten_stack_linear(stack, settings); // Change this later if the chemistry necessitates it.
                write_stack(stack, settings);
                stack->n = 1;
                stack->a[0] = bam_dup1(b);
            } stack_insert(stack, b);
        }
        flatten_stack_linear(stack, settings); // Change this later if the chemistry necessitates it.
        write_stack(stack, settings);
        stack->n = 1;
        bam_destroy1(b);
        // Handle any unpaired reads, though there shouldn't be any in real datasets.
        // What this really means is that the correctly collapsed reads just have a naming issue.
        for(auto pair: settings->realign_pairs) fputs(pair.second.c_str(), settings->fqh);
    }

    void bam_rsq_bookends(rsq_aux_t *settings)
    {

        dlib::tmp_stack_t stack = {0};
        dlib::resize_stack(&stack, STACK_START);
        if(!stack.a)
            LOG_EXIT("Failed to start array of bam1_t structs...\n");
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

        while ((c = getopt(argc, argv, "l:f:t:au?h")) >= 0) {
            switch (c) {
            case 'u':
                settings.cmpkey = UNCLIPPED;
                LOG_INFO("Unclipped start position chosen for cmpkey.\n");
                break;
            case 't': settings.mmlim = atoi(optarg); break;
            case 'f': fqname = optarg; break;
            case 'l': wmode[2] = atoi(optarg)%10 + '0';break;
            case '?': case 'h': return rsq_usage(EXIT_SUCCESS);
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
            for(const char *tag: {"SU", "MU"})
                dlib::check_bam_tag_exit(argv[optind], tag);
        for(const char *tag: {"FM", "FA", "PV", "FP", "RV", "LM"})
            dlib::check_bam_tag_exit(argv[optind], tag);
        settings.in = sam_open(argv[optind], "r");
        settings.hdr = sam_hdr_read(settings.in);

        if (settings.hdr == nullptr || settings.hdr->n_targets == 0)
            LOG_EXIT("input SAM does not have header. Abort!\n");

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
