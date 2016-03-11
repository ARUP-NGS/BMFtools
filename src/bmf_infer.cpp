#include "bmf_infer.h"

namespace BMF {

    std::string BamFisherSet::to_fastq() {
        int i;
        std::string seq;
        seq.resize(len);
        std::vector<uint32_t> agrees(len);
        std::vector<uint32_t> full_quals(len); // igamc calculated
        for(i = 0; i < len; ++i) {
            const int argmaxret = arr_max_u32(phred_sums.data(), i); // 0,1,2,3,4/A,C,G,T,N
            agrees[i] = votes[i * 5 + argmaxret];
            full_quals[i] = pvalue_to_phred(igamc_pvalues(n, LOG10_TO_CHI2(phred_sums[i * 5 + argmaxret])));
            // Mask unconfident base calls
            if(full_quals[i] < 2 || (double)agrees[i] / n < MIN_FRAC_AGREED) {
                seq[i] = 'N';
                max_observed_phreds[i] = '#';
            } else {
            	seq[i] = num2nuc(argmaxret);
            	max_observed_phreds[i] += 33;
            }
        }
        const size_t bufsize = (5 * len + 10); // 6 minimum. Let's give 4 extra just in case?
        std::string pvbuf;
        pvbuf.resize(bufsize);
        kstring_t pvks = {0, bufsize, (char *)pvbuf.data()}; // Does not own!
        std::string fabuf;
        fabuf.resize(bufsize);
        kstring_t faks = {0, bufsize, (char *)pvbuf.data()};
        ksprintf(&pvks, "PV:B:I");
        ksprintf(&faks, "FA:B:I");
        for(int i = 0; i < len; ++i) {
            ksprintf(&pvks, ",%u", full_quals[i]);
            ksprintf(&faks, ",%u", agrees[i]);
        }
        pvbuf.resize(pvks.l);
        fabuf.resize(faks.l);
        std::string ret;
        ret.resize(pvks.l + faks.l + get_name().size() + len * 2 + 32);
        // 32 is for "\n+\n" + "\n" + "FP:i:1\tRV:i:0\n" + "\t" + "\t" + "FM:i:[Up to four digits]"
        stringprintf(ret, "@%s %s\t%s\tFM:i:%i\tFP:i:1\tRV:i:0\n%s\n+\n%s\n",
                     get_name().c_str(), pvbuf.c_str(), fabuf.c_str(), n, seq.c_str(), max_observed_phreds.c_str());
        return ret;
    }
    void BamFisherKing::add_to_hash(infer_aux_t *settings) {
        for(auto& set: sets) {
            auto found = settings->realign_pairs.find(set.second.get_name());
            if(found == settings->realign_pairs.end()) {
                settings->realign_pairs.emplace(set.second.get_name(), set.second.to_fastq());
            } else {
                assert(memcmp(found->second.c_str() + 1,
                              set.second.get_name().c_str(),
                              set.second.get_name().size()) == 0);
                if(set.second.get_is_read1()) {
                    fputs(set.second.to_fastq().c_str(), settings->fqh);
                    fputs(found->second.c_str(), settings->fqh);
                } else {
                    fputs(found->second.c_str(), settings->fqh);
                    fputs(set.second.to_fastq().c_str(), settings->fqh);
                }
            }
        }
    }

    /* OUTLINE FOR INFER
     * 1. multidimensional array (in 1-d) of nucleotide counts.
     * 2. multidimensional array (in 1-d) of phred sums.
     * 3. argmax
     * 4. Write as fastq.
     * 5. ???
     * 6. PROFIT
     *
     * */

    static const std::function<int (bam1_t *, bam1_t *)> fns[4] = {&same_stack_pos, &same_stack_pos_se,
                                                                   &same_stack_ucs, &same_stack_ucs_se};

    std::string get_SO(bam_hdr_t *hdr) {
        char *end, *so_start;
        std::string ret;
        if (strncmp(hdr->text, "@HD", 3) != 0) goto NA;
        if ((end = strchr(hdr->text, '\n')) == 0) goto NA;
        *end = '\0';

        if((so_start = strstr(hdr->text, "SO:")) == nullptr) goto NA;
        ret = std::string(so_start + strlen("SO:"));
        *end = '\n';
        return ret;

        NA:
        LOG_WARNING("Sort order not found. Returning N/A.\n");
        return std::string("N/A");
    }

    inline void bam2ffq(bam1_t *b, FILE *fp)
    {
        char *qual, *seqbuf;
        int i;
        uint8_t *seq, *rvdata;
        uint32_t *pv, *fa;
        int8_t t;
        kstring_t ks = {0, 0, nullptr};
        ksprintf(&ks, "@%s PV:B:I", bam_get_qname(b));
        pv = (uint32_t *)dlib::array_tag(b, "PV");
        fa = (uint32_t *)dlib::array_tag(b, "FA");
        for(i = 0; i < b->core.l_qseq; ++i) ksprintf(&ks, ",%u", pv[i]);
        kputs("\tFA:B:I", &ks);
        for(i = 0; i < b->core.l_qseq; ++i) ksprintf(&ks, ",%u", fa[i]);
        ksprintf(&ks, "\tFM:i:%i\tFP:i:%i\tNC:i:%i",
                bam_itag(b, "FM"), bam_itag(b, "FP"), bam_itag(b, "NC"));
        if((rvdata = bam_aux_get(b, "RV")) != nullptr)
            ksprintf(&ks, "\tRV:i:%i", bam_aux2i(rvdata));
        kputc('\n', &ks);
        seq = bam_get_seq(b);
        seqbuf = (char *)malloc(b->core.l_qseq + 1);
        for (i = 0; i < b->core.l_qseq; ++i) seqbuf[i] = seq_nt16_str[bam_seqi(seq, i)];
        seqbuf[i] = '\0';
        if (b->core.flag & BAM_FREVERSE) { // reverse complement
            for(i = 0; i < b->core.l_qseq>>1; ++i) {
                t = seqbuf[b->core.l_qseq - i - 1];
                seqbuf[b->core.l_qseq - i - 1] = nuc_cmpl(seqbuf[i]);
                seqbuf[i] = nuc_cmpl(t);
            }
            if(b->core.l_qseq&1) seqbuf[i] = nuc_cmpl(seqbuf[i]);
        }
        seqbuf[b->core.l_qseq] = '\0';
        assert(strlen(seqbuf) == (uint64_t)b->core.l_qseq);
        kputs(seqbuf, &ks);
        kputs("\n+\n", &ks);
        qual = (char *)bam_get_qual(b);
        for(i = 0; i < b->core.l_qseq; ++i) seqbuf[i] = 33 + qual[i];
        if (b->core.flag & BAM_FREVERSE) { // reverse
            for (i = 0; i < b->core.l_qseq>>1; ++i) {
                t = seqbuf[b->core.l_qseq - 1 - i];
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

    static inline void update_bam1(bam1_t *p, bam1_t *b)
    {
        uint8_t *bdata, *pdata;
        int n_changed;
        bdata = bam_aux_get(b, "FM");
        pdata = bam_aux_get(p, "FM");
        if(UNLIKELY(!bdata || !pdata)) {
            fprintf(stderr, "Required FM tag not found. Abort mission!\n");
            exit(EXIT_FAILURE);
        }
        int bFM = bam_aux2i(bdata);
        int pFM = bam_aux2i(pdata);
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
        n_changed = dlib::int_tag_zero(pdata) + dlib::int_tag_zero(bdata);
        if(pdata) bam_aux_del(p, pdata);

        uint32_t *const bPV = (uint32_t *)dlib::array_tag(b, "PV"); // Length of this should be b->l_qseq
        uint32_t *const pPV = (uint32_t *)dlib::array_tag(p, "PV"); // Length of this should be b->l_qseq
        uint32_t *const bFA = (uint32_t *)dlib::array_tag(b, "FA"); // Length of this should be b->l_qseq
        uint32_t *const pFA = (uint32_t *)dlib::array_tag(p, "FA"); // Length of this should be b->l_qseq
        uint8_t *const bSeq = (uint8_t *)bam_get_seq(b);
        uint8_t *const pSeq = (uint8_t *)bam_get_seq(p);
        uint8_t *const bQual = (uint8_t *)bam_get_qual(b);
        uint8_t *const pQual = (uint8_t *)bam_get_qual(p);
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


    void write_stack(dlib::tmp_stack_t *stack, infer_aux_t *settings)
    {
        for(unsigned i = 0; i < stack->n; ++i) {
            if(stack->a[i]) {
                uint8_t *data;
                std::string qname;
                if((data = bam_aux_get(stack->a[i], "NC")) != nullptr) {
                    qname = bam_get_qname(stack->a[i]);
                    if(settings->realign_pairs.find(qname) == settings->realign_pairs.end()) {
                        settings->realign_pairs[qname] = dlib::bam2cppstr(stack->a[i]);
                    } else {
                        // Make sure the read names/barcodes match.
                        assert(memcmp(settings->realign_pairs[qname].c_str() + 1, qname.c_str(), qname.size()) == 0);
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
                } else sam_write1(settings->out, settings->hdr, stack->a[i]);
                bam_destroy1(stack->a[i]), stack->a[i] = nullptr;
            }
        }
    }

    static inline void flatten_stack_linear(dlib::tmp_stack_t *stack, infer_aux_t *settings)
    {
        // Sort by read names to make sure that any progressive rescuing ends at the same name.
        // return a ? (b ? (int)(strcmp(bam_get_qname(a), bam_get_qname(b)) < 0): 0): (b ? 1: 0);
        // Returns 0 if comparing two nulls, and returns true that a nullptr lt a valued name
        // Compares strings otherwise.
        BamFisherKing king(stack);
        king.add_to_hash(settings);
    }

    static const char *sorted_order_strings[2] = {"positional_rescue", "unclipped_rescue"};

    void infer_core(infer_aux_t *settings, dlib::tmp_stack_t *stack)
    {
        // This selects the proper function to use for deciding if reads belong in the same stack.
        // It chooses the single-end or paired-end based on is_se and the bmf or pos based on cmpkey.
        std::function<int (bam1_t *, bam1_t *)> fn = fns[settings->is_se | (settings->cmpkey<<1)];
        if(strcmp(get_SO(settings->hdr).c_str(), sorted_order_strings[settings->cmpkey])) {
            LOG_EXIT("Sort order (%s) is not expected %s for rescue mode. Abort!\n",
                     get_SO(settings->hdr).c_str(), sorted_order_strings[settings->cmpkey]);
        }
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
            } else {
                // Keep adding bam records.
                stack_insert(stack, b);
            }
        }
        flatten_stack_linear(stack, settings); // Change this later if the chemistry necessitates it.
        write_stack(stack, settings);
        stack->n = 1;
        bam_destroy1(b);
        // Handle any unpaired reads, though there shouldn't be any in real datasets.
        if(settings->realign_pairs.size()) {
            LOG_WARNING("There shoudn't be any orphaned reads left in real datasets, but there are %lu. Something may be wrong....\n", settings->realign_pairs.size());
            for(auto& pair: settings->realign_pairs) {
                fprintf(settings->fqh, pair.second.c_str());
                settings->realign_pairs.erase(pair.first);
            }
        }
    }

    int infer_usage(int retcode)
    {
        fprintf(stderr,
                        "Positional rescue. \n"
                        "Reads with the same start position (or unclipped start, if -u is set) are compared.\n"
                        "If their barcodes are sufficiently similar, they are treated as having originated"
                        "from the same original template molecule.\n"
                        "Usage:  bmftools infer <input.srt.bam> <output.bam>\n\n"
                        "Flags:\n"
                        "-f      Path for the fastq for reads that need to be realigned. REQUIRED.\n"
                        "-t      Mismatch limit. Default: 2\n"
                        "-l      Set bam compression level. Valid: 0-9. (0 == uncompresed)\n"
                        "-u      Flag to use unclipped start positions instead of pos/mpos for identifying potential duplicates.\n"
                        "Note: -u requires pre-processing with bmftools mark_unclipped.\n"
                );
        return retcode;
    }


    int infer_main(int argc, char *argv[])
    {
        int c;
        char wmode[4] = "wb";

        infer_aux_t settings = {0};
        settings.mmlim = 2;
        settings.cmpkey = POSITION;

        char *fqname = nullptr;

        if(argc < 3) return infer_usage(EXIT_FAILURE);

        while ((c = getopt(argc, argv, "l:f:t:au?h")) >= 0) {
            switch (c) {
            case 'u':
                settings.cmpkey = UNCLIPPED;
                LOG_INFO("Unclipped start position chosen for cmpkey.\n");
                break;
            case 't': settings.mmlim = atoi(optarg); break;
            case 'f': fqname = optarg; break;
            case 'l': wmode[2] = atoi(optarg)%10 + '0';break;
            case '?': case 'h': return infer_usage(EXIT_SUCCESS);
            }
        }
        if (optind + 2 > argc)
            return infer_usage(EXIT_FAILURE);

        if(!fqname) {
            fprintf(stderr, "Fastq path for rescued reads required. Abort!\n");
            return infer_usage(EXIT_FAILURE);
        }

        settings.fqh = fopen(fqname, "w");

        if(!settings.fqh)
            LOG_EXIT("Failed to open output fastq for writing. Abort!\n");

        if(settings.cmpkey == UNCLIPPED)
            for(const char *tag: {"SU", "MU"})
                dlib::check_bam_tag_exit(argv[optind], tag);
        dlib::check_bam_tag_exit(argv[optind], "LM");
        settings.in = sam_open(argv[optind], "r");
        settings.hdr = sam_hdr_read(settings.in);

        if (settings.hdr == nullptr || settings.hdr->n_targets == 0)
            LOG_EXIT("input SAM does not have header. Abort!\n");

        settings.out = sam_open(argv[optind+1], wmode);
        if (settings.in == 0 || settings.out == 0)
            LOG_EXIT("fail to read/write input files\n");
        sam_hdr_write(settings.out, settings.hdr);

        dlib::tmp_stack_t stack = {0};
        dlib::resize_stack(&stack, STACK_START);
        if(!stack.a)
            LOG_EXIT("Failed to start array of bam1_t structs...\n");

        infer_core(&settings, &stack);

        free(stack.a);
        bam_hdr_destroy(settings.hdr);
        sam_close(settings.in); sam_close(settings.out);
        if(settings.fqh) fclose(settings.fqh);
        LOG_INFO("Successfully completed bmftools infer.\n");
        return EXIT_SUCCESS;
    }
}
