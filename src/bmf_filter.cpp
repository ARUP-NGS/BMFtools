#include <getopt.h>
#include <functional>
#include "dlib/bam_util.h"

namespace bmf {

    int usage(char **argv, int retcode=EXIT_FAILURE) {
        fprintf(stderr,
                        "Filters a bam by a set of given parameters.\n"
                        "Usage: bmftools filter <-l output_compression_level> in.bam out.bam\n"
                        "Use - for stdin or stdout.\n"
                        "Flags:\n"
                        "-m\t\tFail reads with mapping quality < parameter.\n"
                        "-a\t\tFail read pairs where both reads' aligned fraction < parameter.\n"
                        "-F\t\tSkip all reads with any bits in parameter set.\n"
                        "-f\t\tSkip reads not sharing all bits in parameter set.\n"
                        "-b\t\tPath to bed file with which to filter.\n"
                        "-P\t\tNumber of bases around bed file regions to pad.\n"
                        "-s\t\tMinimum family size for inclusion.\n"
                        "-r\t\tIf set, writes failed reads to this file.\n"
                        "-v\t\tInvert pass/fail, analogous to grep.\n"
                );
        return retcode;
    }

    struct opts {
        uint32_t minFM:14;
        uint32_t minmq:8;
        uint32_t v:1;
        uint32_t is_se:1;
        uint32_t skip_flag:16;
        uint32_t require_flag:16;
        float minAF;
        khash_t(bed) *bed;
    };

/* If FM tag absent, it's treated as if it were 1.
 * Fail reads with FM < minFM, MQ < minmq, a flag with any skip bits set,
 * a flag without all required bits set, and, if a bed file is provided,
 * reads outside of the bed region.
*/
    static inline int test_core(bam1_t *b, opts *options) {
        uint8_t *data;
        if((((data = bam_aux_get(b, "FM")) != nullptr) ? bam_aux2i(data): 1) >= (int)(options->minFM))
            if(b->core.qual >= options->minmq)
                if((b->core.flag & options->skip_flag) == 0)
                    if((b->core.flag & options->require_flag) == options->require_flag)
                        if(options->bed ? dlib::bed_test(b, options->bed):1)
                            if(((data = bam_aux_get(b, "MF")) == nullptr ? 1: bam_aux2i(data) >= options->minAF)
                               || dlib::bam_frac_align(b) >= options->minAF)
                                return 1;
        return 0;
    }

    /* If FM tag absent, it's treated as if it were 1.
     * Fail reads with FM < minFM, MQ < minmq, a flag with any skip bits set,
     * a flag without all required bits set, and, if a bed file is provided,
     * reads outside of the bed region.
    */
        static inline int test_core_se(bam1_t *b, opts *options) {
            uint8_t *data;
            return ((((data = bam_aux_get(b, "FM")) != nullptr) ? bam_aux2i(data): 1) >= (int)(options->minFM)) &&
                    b->core.qual >= options->minmq &&
                    ((b->core.flag & options->skip_flag) == 0) &&
                    (b->core.flag & options->require_flag) == options->require_flag &&
                    (options->bed ? dlib::bed_test(b, options->bed):1) &&
                    dlib::bam_frac_align(b) >= options->minAF;
        }

    /*
     * Return a 0 status to pass, a 1 to fail, in order to match for_each.
     */
    int bam_test(bam1_t *b, void *options) {
        if(((opts *)options)->is_se)
            return ((opts *)options)->v ? test_core_se(b, (opts *)options)
                                        : !test_core_se(b, (opts *)options);
        return ((opts *)options)->v ? test_core(b, (opts *)options)
                                    : !test_core(b, (opts *)options);
    }

    int filter_split_core(dlib::BamHandle& in, dlib::BamHandle& out, dlib::BamHandle& refused, opts *param)
    {
        uint64_t count = 0;
        while(in.next() >= 0) {
            if(++count % 1000000 == 0) LOG_INFO("%lu records processed.\n", count);
            if(bam_test(in.rec, (void *)param)) {
                refused.write(in.rec);
            } else {
                out.write(in.rec);
            }
        }
        return EXIT_SUCCESS;
    }

    int filter_main(int argc, char *argv[]) {
        if(argc < 3)
            return usage(argv);
        if(strcmp(argv[1], "--help") == 0)
            return usage(argv, EXIT_SUCCESS);
        int c;
        char out_mode[4]{"wb"};
        opts param = {0};
        char *bedpath = nullptr;
        int padding = DEFAULT_PADDING;
        std::string refused_path("");
        while((c = getopt(argc, argv, "s:a:r:P:b:m:F:f:l:hAv?")) > -1) {
            switch(c) {
            case 'a': param.minAF = atof(optarg); break;
            case 'P': padding = atoi(optarg); break;
            case 'b': bedpath = optarg; break;
            case 'm': param.minmq = strtoul(optarg, nullptr, 0); break;
            case 's': param.minFM = strtoul(optarg, nullptr, 0); break;
            case 'F': param.skip_flag = strtoul(optarg, nullptr, 0); break;
            case 'f': param.require_flag = strtoul(optarg, nullptr, 0); break;
            case 'v': param.v = 1; break;
            case 'r': refused_path = optarg; break;
            case 'l': out_mode[2] = atoi(optarg) % 10 + '0'; break;
            case '?': case 'h': return usage(argv, EXIT_SUCCESS);
            }
        }
        dlib::check_bam_tag_exit(argv[optind], "FM");
        if(argc - 2 != optind) {
            LOG_EXIT("Required: precisely two positional arguments (in bam, out bam).\n");
        }
        dlib::BamHandle in(argv[optind]);
        param.bed = bedpath ? dlib::parse_bed_hash(bedpath, in.header, padding)
                            : nullptr;
        dlib::add_pg_line(in.header, argc, argv, "bmftools filter", BMF_VERSION,
                "bmftools", "Filters or splits a bam by a set of criteria.");
        if(param.minAF > 0 && param.is_se == 0)
            dlib::check_bam_tag_exit(argv[optind], "MF");
        dlib::BamHandle out(argv[optind + 1], in.header, out_mode);
        // Core
        int ret = -1;
        if(refused_path.size()) { // refused path is set.
            LOG_DEBUG("Writing passing records to %s, failing to %s.\n", out.fp->fn, refused_path.c_str());
            dlib::BamHandle refused(refused_path.c_str(), in.header, out_mode);
            ret = filter_split_core(in, out, refused, &param);
        } else ret = in.for_each(bam_test, out, (void *)&param);
        // Clean up.
        if(param.bed) dlib::bed_destroy_hash((void *)param.bed);
        LOG_INFO("Successfully completed bmftools filter!\n");
        return ret;
    }

}
