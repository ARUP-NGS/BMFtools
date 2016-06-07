#include "bmf_stack.h"

#include <getopt.h>
#include <algorithm>

namespace bmf {

    void stack_usage(int retcode)
    {
        fprintf(stderr,
                        "Builds a stack summary for base calls and performs simple filters"
                        " to produce a maximally permissive 'variant caller'.\n"
                        "Usage:\nbmftools stack <opts> <tumor.srt.indexed.bam> <normal.srt.indexed.bam>\n"
                        "Optional arguments:\n"
                        "-R, --ref\tPath to fasta reference. REQUIRED.\n"
                        "-o, --outpath\tPath to output file. Defaults to stdout.\n"
                        "-b, --bedpath\tPath to bed file to only validate variants in said region. REQUIRED.\n"
                        "-c, --min-count\tMinimum number of observations for a given allele passing filters to pass variant. Default: 1.\n"
                        "-s, --min-family-size\tMinimum number of reads in a family to include a that collapsed observation\n"
                        "-D, --min-duplex\tMinimum number of duplex reads supporting a variant to pass it. Default: 0.\n"
                        "-O, --min-overlap\tMinimum number of concordant overlapping read-pairs supporting a variant to pass it. Default: 0.\n"
                        "-P, --skip-improper\tSkip reads not marked as being in a proper pair.\n"
                        "-2, --skip-secondary\tSkip secondary alignments.\n"
                        "-S, --skip-supplementary\tSkip supplementary alignments.\n"
                        "-q, --skip-qcfail\tSkip reads marked as QC fail.\n"
                        "-r, --skip-duplicates\tSkip reads marked as being PCR or optical duplicates.\n"
                        "-F, --skip-recommended\tSkip reads marked as QC fail, duplicate, or secondary. Equivalent to -P -r -2 -q\n"
                        "-f, --min-fraction-agreed\tMinimum fraction of reads in a family agreed on a base call\n"
                        "-v, --min-phred-quality\tMinimum calculated p-value on a base call in phred space\n"
                        "-p, --padding\tNumber of bases outside of bed region to pad.\n"
                        "-a, --min-family-agreed\tMinimum number of reads in a family agreed on a base call\n"
                        "-m, --min-mapping-quality\tMinimum mapping quality for reads for inclusion\n"
                        "-B, --emit-bcf-format\tEmit bcf-formatted output. (Defaults to vcf).\n"
                );
        exit(retcode);
    }


    static int read_bam(dlib::BamHandle *data, bam1_t *b)
    {
        int ret;
        uint8_t *tmp;
        for(;;)
        {
            if(!data->iter) LOG_EXIT("Need to access bam with index.\n");
            ret = sam_itr_next(data->fp, data->iter, b);
            if ( ret<0 ) break;
            if((tmp = bam_aux_get(b, "FP")) != nullptr && bam_aux2i(tmp))
                if((b->core.flag & BAM_FUNMAP) == 0)
                    break;
        }
        return ret;
    }

    void process_matched_pileups(bmf::stack_aux_t *aux, bcf1_t *ret,
                            const int& tn_plp, const int& tpos, const int& ttid,
                            const int& nn_plp, const int& npos, const int& ntid) {
        // Build overlap hash
        std::unordered_map<std::string, bmf::UniqueObservation> tobs, nobs;
        std::unordered_map<std::string, bmf::UniqueObservation>::iterator found;
        int flag_failed[2]{0};
        int af_failed[2]{0};
        int mq_failed[2]{0};
        int improper_count[2]{0};
        int olap_count[2]{0};
        std::string qname;
        for(int i = 0; i < tn_plp; ++i) {
            if(aux->tumor.pileups[i].is_del || aux->tumor.pileups[i].is_refskip) continue;
            if(aux->conf.skip_flag & aux->tumor.pileups[i].b->core.flag) {
                ++flag_failed[0];
                continue;
            }
            if((aux->tumor.pileups[i].b->core.flag & BAM_FPROPER_PAIR) == 0) {
                ++improper_count[0];
                if(aux->conf.skip_improper) continue;
            }
            if(dlib::bam_frac_align(aux->tumor.pileups[i].b) < aux->conf.minAF) {
                ++af_failed[0]; continue;
            }
            // Add in

            qname = bam_get_qname(aux->tumor.pileups[i].b);
            if((found = tobs.find(qname)) == tobs.end()) {
                //LOG_DEBUG("Put in entry at index %i with tn_plp as %i\n", i, tn_plp);
                tobs.emplace(qname, aux->tumor.pileups[i]);
            } else {
                ++olap_count[0];
                //LOG_DEBUG("Added other in pair with qname %s.\n", qname.c_str());
                found->second.add_obs(aux->tumor.pileups[i]);
            }
        }
        for(auto& pair: tobs) {
            if(pair.second.get_max_mq() < aux->conf.minmq)
                ++mq_failed[0], pair.second.set_pass(0);
        }
        for(int i = 0; i < nn_plp; ++i) {
            if(aux->normal.pileups[i].is_del || aux->normal.pileups[i].is_refskip) continue;
            if(aux->conf.skip_flag & aux->normal.pileups[i].b->core.flag) {
                ++flag_failed[1];
                continue;
            }
            if((aux->normal.pileups[i].b->core.flag & BAM_FPROPER_PAIR) == 0) {
                ++improper_count[1];
                if(aux->conf.skip_improper) continue;
            }
            if(dlib::bam_frac_align(aux->normal.pileups[i].b) < aux->conf.minAF) {
                ++af_failed[1]; continue;
            }
            qname = bam_get_qname(aux->normal.pileups[i].b);
            if((found = nobs.find(qname)) == nobs.end())
                nobs.emplace(qname, aux->normal.pileups[i]);
            else ++olap_count[1], found->second.add_obs(aux->normal.pileups[i]);
        }
        for(auto& pair: nobs) {
            if(pair.second.get_max_mq() < aux->conf.minmq)
                ++mq_failed[1], pair.second.set_pass(0);
        }
        //LOG_DEBUG("Making PairVCFPos.\n");
        // Build vcfline struct
        bmf::PairVCFPos vcfline(tobs, nobs, ttid, tpos);
        vcfline.to_bcf(ret, aux, ttid, tpos);
        bcf_update_format_int32(aux->vcf.vh, ret, "MQ_FAILED", (void *)mq_failed, COUNT_OF(mq_failed) * 2);
        bcf_update_format_int32(aux->vcf.vh, ret, "AF_FAILED", (void *)af_failed, COUNT_OF(af_failed) * 2);
        bcf_update_format_int32(aux->vcf.vh, ret, "OVERLAP", (void *)olap_count, COUNT_OF(olap_count) * 2);
        //LOG_INFO("Ret for writing vcf to file: %i.\n", aux->vcf.write(ret));
        aux->vcf.write(ret);
        bcf_clear(ret);
    }

    /*
     * Needs a rewrite after the T/N pair rewrite!
     */
    void process_pileup(bcf1_t *ret, const bam_pileup1_t *plp, int n_plp, int pos, int tid, bmf::stack_aux_t *aux) {
        // Build overlap hash
        std::unordered_map<std::string, bmf::UniqueObservation> obs;
        std::unordered_map<std::string, bmf::UniqueObservation>::iterator found;
        int flag_failed(0);
        int af_failed(0);
        int mq_failed(0);
        int improper_count(0);
        int olap_count(0);
        std::string qname;
        for(int i = 0; i < n_plp; ++i) {
            if(aux->tumor.pileups[i].is_del || aux->tumor.pileups[i].is_refskip) continue;
            if(aux->conf.skip_flag & aux->tumor.pileups[i].b->core.flag) {
                ++flag_failed;
                continue;
            }
            if((aux->tumor.pileups[i].b->core.flag & BAM_FPROPER_PAIR) == 0) {
                ++improper_count;
                if(aux->conf.skip_improper) continue;
            }
            if(dlib::bam_frac_align(aux->tumor.pileups[i].b) < aux->conf.minAF) {
                ++af_failed; continue;
            }
            qname = bam_get_qname(aux->tumor.pileups[i].b);
            if((found = obs.find(qname)) == obs.end()) {
                //LOG_DEBUG("Put in entry at index %i with tn_plp as %i\n", i, tn_plp);
                obs.emplace(qname, aux->tumor.pileups[i]);
            } else {
                ++olap_count;
                //LOG_DEBUG("Added other in pair with qname %s.\n", qname.c_str());
                found->second.add_obs(aux->tumor.pileups[i]);
            }
        }
        for(auto& pair: obs)
            if(pair.second.get_max_mq() < aux->conf.minmq)
                ++mq_failed, pair.second.set_pass(0);
        //LOG_DEBUG("Making PairVCFPos.\n");
        // Build vcfline struct
        bmf::SampleVCFPos vcfline(obs, tid, pos);
        vcfline.to_bcf(ret, aux, aux->get_ref_base(tid, pos));
        bcf_update_info_int32(aux->vcf.vh, ret, "MQ_FAILED", (void *)&mq_failed, 1);
        bcf_update_info_int32(aux->vcf.vh, ret, "AF_FAILED", (void *)&af_failed, 1);
        bcf_update_info_int32(aux->vcf.vh, ret, "IMPROPER", (void *)&improper_count, 1);
        bcf_update_format_int32(aux->vcf.vh, ret, "OVERLAP", (void *)&olap_count, 1);
        aux->vcf.write(ret);
        bcf_clear(ret);
    }

    int stack_core_single(bmf::stack_aux_t *aux)
    {
        if(!aux->tumor.idx)
            LOG_EXIT("Could not load bam index. Abort!\n");
        aux->tumor.plp = bam_plp_init((bam_plp_auto_f)read_bam, (void *)&aux->tumor);
        LOG_DEBUG("Max depth: %i.\n", aux->conf.max_depth);
        bam_plp_set_maxcnt(aux->tumor.plp, aux->conf.max_depth);
        LOG_DEBUG("Making sorted keys.\n");
        std::vector<khiter_t> sorted_keys(dlib::make_sorted_keys(aux->bed));
        int tid, pos, n_plp;
        bcf1_t *v(bcf_init1());
        for(unsigned k = 0; k < sorted_keys.size(); ++k) {
            const khiter_t key = sorted_keys[k];
            LOG_DEBUG("Now iterating through tid %i.\n", kh_key(aux->bed, key));
            const size_t n = kh_val(aux->bed, key).n;
            for(uint64_t i = 0; i < n; ++i) {
                const int start(get_start(kh_val(aux->bed, key).intervals[i]));
                const int stop(get_stop(kh_val(aux->bed, key).intervals[i]));
                const int bamtid((int)kh_key(aux->bed, key));
                if(aux->single_region_itr(bamtid, start, stop, n_plp, pos, tid))
                    continue;  // Could not load reads in one of the two bams.
                process_pileup(v, aux->tumor.pileups, n_plp, pos, tid, aux);
                while(aux->next_single_pileup(&tid, &pos, &n_plp, stop))
                    process_pileup(v, aux->tumor.pileups, n_plp, pos, tid, aux);
            }
        }
        bcf_destroy(v);
        return 0;
    }

    int stack_core(bmf::stack_aux_t *aux)
    {
        if(!aux->tumor.idx || !aux->normal.idx)
            LOG_EXIT("Could not load bam indices. Abort!\n");
        aux->tumor.plp = bam_plp_init((bam_plp_auto_f)read_bam, (void *)&aux->tumor);
        aux->normal.plp = bam_plp_init((bam_plp_auto_f)read_bam, (void *)&aux->normal);
        LOG_DEBUG("Max depth: %i.\n", aux->conf.max_depth);
        bam_plp_set_maxcnt(aux->tumor.plp, aux->conf.max_depth);
        bam_plp_set_maxcnt(aux->normal.plp, aux->conf.max_depth);
        LOG_DEBUG("Making sorted keys.\n");
        std::vector<khiter_t> sorted_keys(dlib::make_sorted_keys(aux->bed));
        int ttid, tpos, tn_plp, ntid, npos, nn_plp;
        bcf1_t *v = bcf_init1();
        for(unsigned k = 0; k < sorted_keys.size(); ++k) {
            const khiter_t key = sorted_keys[k];
            LOG_DEBUG("Now iterating through tid %i.\n", kh_key(aux->bed, key));
            const size_t n = kh_val(aux->bed, key).n;
            for(uint64_t i = 0; i < n; ++i) {
                const int start = get_start(kh_val(aux->bed, key).intervals[i]);
                const int stop = get_stop(kh_val(aux->bed, key).intervals[i]);
                const int bamtid = (int)kh_key(aux->bed, key);
                if(aux->pair_region_itr(bamtid, start, stop, tn_plp, tpos, ttid, nn_plp, npos, ntid))
                    continue;  // Could not load reads in one of the two bams.
                process_matched_pileups(aux, v, tn_plp, tpos, ttid, nn_plp, npos, ntid);
                while(aux->next_paired_pileup(&ttid, &tpos, &tn_plp, &ntid, &npos, &nn_plp, stop))
                    process_matched_pileups(aux, v, tn_plp, tpos, ttid, nn_plp, npos, ntid);
            }
        }
        bcf_destroy(v);
        return 0;
    }

    int stack_main(int argc, char *argv[]) {
        int c;
        unsigned padding = (unsigned)-1;
        if(argc < 2) stack_usage(EXIT_FAILURE);
        char *outvcf = nullptr, *refpath = nullptr;
        char *bedpath = nullptr;
        struct bmf::stack_conf_t conf = {0};
        const struct option lopts[] = {
            {"skip-secondary", no_argument, nullptr, '2'},
            {"min-family-agreed", required_argument, nullptr, 'a'},
            {"bedpath", required_argument, nullptr, 'b'},
            {"emit-bcf", no_argument, nullptr, 'B'},
            {"min-count", required_argument, nullptr, 'c'},
            {"max-depth", required_argument, nullptr, 'd'},
            {"min-duplex", required_argument, nullptr, 'D'},
            {"min-fraction-agreed", required_argument, nullptr, 'f'},
            {"min-mapping-quality", required_argument, nullptr, 'm'},
            {"min-overlap", required_argument, nullptr, 'O'},
            {"out-vcf", required_argument, nullptr, 'o'},
            {"skip-improper", no_argument, nullptr, 'P'},
            {"padding", required_argument, nullptr, 'p'},
            {"skip-qcfail", no_argument, nullptr, 'q'},
            {"skip-duplicates", no_argument, nullptr, 'r'},
            {"skip-recommended", no_argument, nullptr, 'F'},
            {"ref", required_argument, nullptr, 'R'},
            {"min-family-size", required_argument, nullptr, 's'},
            {"skip-supplementary", no_argument, nullptr, 'S'},
            {"min-phred-quality", required_argument, nullptr, 'v'},
            {0, 0, 0, 0}
        };
        while ((c = getopt_long(argc, argv, "R:D:q:r:2:S:d:a:s:m:p:f:b:v:o:O:c:BP?hVF", lopts, nullptr)) >= 0) {
            switch (c) {
                case '2': conf.skip_flag |= BAM_FSECONDARY; break;
                case 'a': conf.minFA = atoi(optarg); break;
                case 'b': bedpath = optarg; break;
                case 'B': conf.output_bcf = 1; break;
                case 'c': conf.min_count = atoi(optarg); break;
                case 'd': conf.max_depth = atoi(optarg); break;
                case 'D': conf.min_duplex = atoi(optarg); break;
                case 'F':
                          conf.skip_flag |= (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP);
                          conf.skip_improper = 1;
                          break;
                case 'f': conf.min_fr = (float)atof(optarg); break;
                case 'm': conf.minmq = atoi(optarg); break;
                case 'O': conf.min_overlap = atoi(optarg); break;
                case 'o': outvcf = optarg; break;
                case 'P': conf.skip_improper = 1; break;
                case 'p': padding = atoi(optarg); break;
                case 'q': conf.skip_flag |= BAM_FQCFAIL; break;
                case 'r': conf.skip_flag |= BAM_FDUP; break;
                case 'R': refpath = optarg; break;
                case 's': conf.minFM = atoi(optarg); break;
                case 'S': conf.skip_flag |= BAM_FSUPPLEMENTARY; break;
                case 'v': conf.minPV = atoi(optarg); break;
                case 'h': case '?': stack_usage(EXIT_SUCCESS);
            }
        }
        /*
        if(!conf.max_depth) {
            conf.max_depth = DEFAULT_MAX_DEPTH;
        }
        */
        if(optind == argc - 1) {
            LOG_INFO("One bam provided. Running in se mode.\n");
        } else if(optind > argc - 1) {
            LOG_EXIT("Insufficient arguments. Input bam required!\n");
        }
        if(padding < 0) {
            LOG_WARNING("Padding not set. Using default %i.\n", DEFAULT_PADDING);
            padding = DEFAULT_PADDING;
        }
        if(!refpath) {
            LOG_EXIT("refpath required. Abort!\n");
        }
        if(!outvcf) outvcf = (char *)"-";
        bcf_hdr_t *vh = bcf_hdr_init(conf.output_bcf ? "wb": "w");
        add_stack_lines(vh);
        // Add samples
        int tmp;
        if((tmp = bcf_hdr_add_sample(vh, "Tumor")))
            LOG_EXIT("Could not add name %s. Code: %i.\n", "Tumor", tmp);
        if((tmp = bcf_hdr_add_sample(vh, "Normal")))
            LOG_EXIT("Could not add name %s. Code: %i.\n", "Normal", tmp);
        // Add header lines
        bcf_hdr_add_sample(vh, nullptr);
        bcf_hdr_nsamples(vh) = 2;
        // Add command line call
        kstring_t tmpstr = {0};
        ksprintf(&tmpstr, "##cmdline=");
        kputs("bmftools ", &tmpstr);
        for(int i = 0; i < argc; ++i) ksprintf(&tmpstr, "%s ", argv[i]);
        bcf_hdr_append(vh, tmpstr.s);
        tmpstr.l = 0;
        bcf_hdr_printf(vh, "##bed_filename=\"%s\"", bedpath ? bedpath: "FullGenomeAnalysis");
        samFile *tmpfp = sam_open(argv[optind], "r");
        bam_hdr_t *hdr = sam_hdr_read(tmpfp);
        free(tmpstr.s);
        std::string timestring("", 16uL);
        dlib::string_fmt_time(timestring);
        bcf_hdr_printf(vh, "##StartTime=\"%s\"", timestring.c_str());
        dlib::bcf_add_bam_contigs(vh, hdr);
        bmf::stack_aux_t aux(argv[optind], argc - 1 == optind ? nullptr: argv[optind + 1], outvcf, vh, conf);
        bcf_hdr_destroy(vh);
        bam_hdr_destroy(hdr);
        if((aux.fai = fai_load(refpath)) == nullptr) LOG_EXIT("failed to open fai. Abort!\n");
        LOG_DEBUG("Bedpath: %s.\n", bedpath);
        if(bedpath == nullptr) LOG_EXIT("Bed path for analysis required.\n");
        aux.bed = dlib::parse_bed_hash(bedpath, aux.tumor.header, padding);
        if(!aux.bed) LOG_EXIT("Could not open bedfile.\n");
        // Check for required tags.
        for(auto tag: {"FM", "FA", "PV", "FP"})
            dlib::check_bam_tag_exit(aux.tumor.fp->fn, tag);
        int ret = (argc - 1 == optind) ? stack_core_single(&aux): stack_core(&aux);
        LOG_INFO("Successfully complete bmftools stack!\n");
        return ret;
    }

} /* namespace bmf */
