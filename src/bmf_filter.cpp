#include "dlib/bam_util.h"
#include "dlib/bed_util.h"
#include <getopt.h>
#include <functional>

namespace BMF {

    int usage(char **argv, int retcode=EXIT_FAILURE) {
        fprintf(stderr,
                        "Filters a bam by a set of given parameters.\n"
                        "Usage: bmftools filter <-l output_compression_level> in.bam out.bam\n"
                        "Use - for stdin or stdout.\n"
                        "Flags:\n"
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
        uint32_t minMQ:8;
        uint32_t v:1;
        uint32_t skip_flag:16;
        uint32_t require_flag:16;
        float minAF;
        khash_t(bed) *bed;
    };

/* If FM tag absent, it's treated as if it were 1.
 * Fail reads with FM < minFM, MQ < minMQ, a flag with any skip bits set,
 * a flag without all required bits set, and, if a bed file is provided,
 * reads outside of the bed region.
*/
    static inline int test_core(bam1_t *b, opts *options) {
        uint8_t *data;
        return ((((data = bam_aux_get(b, "FM")) != nullptr) ? bam_aux2i(data): 1) >= (int)(options->minFM)) &&
                b->core.qual >= options->minMQ &&
                ((b->core.flag & options->skip_flag) == 0) &&
                (b->core.flag & options->require_flag) == options->require_flag &&
                (options->bed ? dlib::bed_test(b, options->bed):1);
    }

    /*
     * Return a 0 status to pass, a 1 to fail, in order to match for_each.
     */
    int bam_test(bam1_t *b, void *options) {
        return ((opts *)options)->v ? test_core(b, (opts *)options)
                                    : !test_core(b, (opts *)options);
    }

    int filter_split_core(dlib::BamHandle& in, dlib::BamHandle& out, dlib::BamHandle& refused, opts *param)
    {
        while(in.next() >= 0) {
            if(bam_test(in.rec, (void *)param)) {
                LOG_DEBUG("Writing to output file %s.\n", refused.fp->fn);
                refused.write(in.rec);
            } else {
                LOG_DEBUG("Writing to output file %s.\n", out.fp->fn);
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
        char out_mode[4] = "wb";
        opts param = {0};
        char *bedpath = nullptr;
        int padding = DEFAULT_PADDING;
        std::string refused_path("");
        while((c = getopt(argc, argv, "s:a:r:P:b:m:F:f:l:hv?")) > -1) {
            switch(c) {
            case 'a': param.minAF = atof(optarg); break;
            case 'P': padding = atoi(optarg); break;
            case 'b': bedpath = optarg; break;
            case 'm': param.minMQ = strtoul(optarg, nullptr, 0); break;
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
        if(bedpath) {
            param.bed = dlib::parse_bed_hash(bedpath, in.header, padding);
        } else {
            LOG_DEBUG("param.bed pointer: %p.\n", (void *)param.bed);
        }
        dlib::BamHandle out(argv[optind + 1], in.header, out_mode);
        // Core
        int ret = -1;
        if(refused_path.size()) { // refused path is set.
            LOG_DEBUG("Splitting. Refused go to %s, pass to %s.\n", refused_path.c_str(), out.fp->fn);
            dlib::BamHandle refused(refused_path.c_str(), in.header, out_mode);
            ret = filter_split_core(in, out, refused, &param);
        } else ret = in.for_each(bam_test, out, (void *)&param);
        // Clean up.
        if(param.bed) dlib::bed_destroy_hash((void *)param.bed);
        LOG_INFO("Successfully completed bmftools filter!\n");
        return ret;
    }

}
