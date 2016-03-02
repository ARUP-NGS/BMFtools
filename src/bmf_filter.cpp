#include "dlib/bam_util.h"
#include "dlib/bed_util.h"
#include <getopt.h>
#include <functional>

int usage(char **argv, int retcode=EXIT_FAILURE) {
    fprintf(stderr, "bmftools %s <-l output_compression_level> in.bam out.bam\n"
                    "Use - for stdin or stdout.\n"
                    "Flags-\n"
                    "-F\t\tSkip all reads with any bits in parameter set.\n"
                    "-f\t\tSkip reads not sharing all bits in parameter set.\n"
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
    float minAF;
    khash_t(bed) *bed;
};

#define TEST(b, data, options) \
        (\
             (\
            (((data = bam_aux_get(b, "FM")) != NULL) ? bam_aux2i(data): 1)\
                                                     >= (int)((opts *)options)->minFM) &&\
             b->core.qual >= ((opts *)options)->minMQ &&\
             ((b->core.flag & ((opts *)options)->skip_flag) == 0) &&\
             ((b->core.flag & ((opts *)options)->require_flag) == (((opts *)options)->require_flag)) &&\
             (((opts *)options)->bed ? bed_test(b, ((opts *)options)->bed)\
                                     : 1)\
        )

#if !NDEBUG
int slow_test(bam1_t *b, uint8_t *data, void *options) {
    data = bam_aux_get(b, "FM");
    if(data) {
        if(bam_aux2i(data) < ((opts *)options)->minFM) {
            //LOG_DEBUG("FM fail (%i).\n", bam_aux2i(data));
            return 0;
        }
    }
    if(b->core.qual < ((opts *)options)->minMQ) {
        //LOG_DEBUG("MQ fail.\n");
        return 0;
    }
    if(b->core.flag & ((opts *)options)->skip_flag) {
        //LOG_DEBUG("Skip flag fail.\n");
        return 0;
    }
    if((b->core.flag & ((opts *)options)->require_flag) != ((opts *)options)->require_flag) {
        //LOG_DEBUG("Require flag fail.\n");
        return 0;
    }
    if(((opts *)options)->bed && !bed_test(b, ((opts *)options)->bed)) {
        //LOG_DEBUG("Bed fail.\n");
        return 0;
    }
    return 1;
}
#endif

int bam_test(bam1_t *b, void *options) {
    uint8_t *data;
#if !NDEBUG
    // Make gcc happy about -Wsequence-point
    int tmp = TEST(b, data, options);
    assert(tmp == slow_test(b, data, options));
#endif
    return ((opts *)options)->v ? !TEST(b, data, options): TEST(b, data, options);
}

#undef TEST

int filter_split_core(dlib::BamHandle& in, dlib::BamHandle& out, dlib::BamHandle& refused, opts *param)
{
    while(in.next() >= 0) (bam_test(in.rec, (void *)param) ? out: refused).write(in.rec);
    /*
    while(in.next() >= 0) {
        if(bam_test(in.rec, (void *)param)) {
            out.write(in.rec);
        } else {
            refused.write(in.rec);
        }
    }
    */
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
    char *bedpath = NULL;
    int padding = DEFAULT_PADDING;
    std::string refused_path("");
    while((c = getopt(argc, argv, "s:a:r:P:b:m:F:f:l:hv?")) > -1) {
        switch(c) {
        case 'a': param.minAF = atof(optarg); break;
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
    } else {
        LOG_DEBUG("param.bed pointer: %p.\n", (void *)param.bed);
    }
    dlib::BamHandle out(argv[optind + 1], in.header, out_mode);
    // Core
    if(refused_path.size()) { // refused path is set.
        LOG_DEBUG("Splitting. Refused go to %s, pass to %s.\n", refused_path.c_str(), out.fp->fn);
        dlib::BamHandle refused(refused_path.c_str(), in.header, out_mode);
        filter_split_core(in, out, refused, &param);
    } else in.for_each(bam_test, out, (void *)&param);
    // Clean up.
    if(param.bed) bed_destroy_hash((void *)param.bed);
    return EXIT_SUCCESS;
}

