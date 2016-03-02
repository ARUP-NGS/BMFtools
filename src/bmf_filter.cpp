#include "dlib/bam_util.h"
#include "dlib/bed_util.h"
#include <getopt.h>
#include <functional>

int usage(char **argv, int retcode=EXIT_FAILURE) {
    fprintf(stderr, "bmftools %s <-l output_compression_level> in.bam out.bam\n"
                    "Use - for stdin or stdout.\n"
                    "Flags-\n"
                    "-F\t\tSkip all reads with any bits in parameter set.\n"
                    "-f\t\tSkip reads sharing no bits in parameter set.\n"
                    "-b\t\tPath to bed file with which to filter.\n"
                    "-P\t\tNumber of bases around bed file regions to pad.\n"
                    "-s\t\tMinimum family size for inclusion.\n"
                    "-r\t\tIf set, writes failed reads to this file.\n"
                    "-v\t\tInvert pass/fail, analogous to grep.\n"
            , argv[0]);
    return retcode;
}

struct opts {
    uint32_t minFM:14;
    uint32_t v:1;
    uint32_t skip_flag:16;
    uint32_t require_flag:16;
    uint32_t minMQ:8;
    khash_t(bed) *bed;
};

#define TEST(b, data, options) \
        (\
             (((data = bam_aux_get(b, "FM")) != NULL) ? bam_aux2i(data) : 1 >= (int)((opts *)options)->minFM) &&\
             b->core.qual >= ((opts *)options)->minMQ &&\
             ((b->core.flag & ((opts *)options)->skip_flag) == 0) &&\
             ((b->core.flag & ((opts *)options)->require_flag)) &&\
             bed_test(b, ((opts *)options)->bed)\
        )

int bam_test(bam1_t *b, void *options) {
    uint8_t *data;
    return ((opts *)options)->v ? !TEST(b, data, options): TEST(b, data, options);
}

#undef TEST

int filter_main(int argc, char *argv[]) {
    if(argc < 3)
        return usage(argv);
    if(strcmp(argv[1], "--help") == 0)
        return usage(argv, EXIT_SUCCESS);
    int c;
    char out_mode[4] = "wb";
    opts param = {0};
    char *bedpath = NULL;
    int padding = DEFAULT_PADDING;
    std::string refused_path("");
    while((c = getopt(argc, argv, "r:P:b:m:F:f:l:hv?")) > -1) {
        switch(c) {
        case 'P': padding = atoi(optarg); break;
        case 'b': bedpath = optarg; break;
        case 'm': param.minMQ = strtoul(optarg, NULL, 0); break;
        case 's': param.minFM = strtoul(optarg, NULL, 0); break;
        case 'F': param.skip_flag = strtoul(optarg, NULL, 0); break;
        case 'f': param.require_flag = strtoul(optarg, NULL, 0); break;
        case 'v': param.v = 1; break;
        case 'r': refused_path = optarg; break;
        case 'l': out_mode[2] = atoi(optarg) % 10 + '0'; break;
        case 'h': case '?': return usage(argv, EXIT_SUCCESS);
        }
    }
    check_bam_tag_exit(argv[optind], "FM");
    if(argc - 2 != optind) {
        LOG_EXIT("Required: precisely two positional arguments (in bam, out bam).\n");
    }
    dlib::BamHandle in(argv[optind]);
    if(bedpath) {
        param.bed = parse_bed_hash(bedpath, in.header, padding);
    }
    dlib::BamHandle out(argv[optind] + 1, in.header, out_mode);
    // Core
    if(!refused_path.size()) {
        in.for_each(bam_test, out, (void *)&param);
    } else {
        dlib::BamHandle refused(refused_path.c_str(), in.header, out_mode);
        while(in.next() >= 0) {
            if(bam_test(in.rec, (void *)&param))
                out.write(in.rec);
            else
                refused.write(in.rec);
        }
    }
    // Clean up.
    if(param.bed) bed_destroy_hash((void *)param.bed);
    return EXIT_SUCCESS;
}

