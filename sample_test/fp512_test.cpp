#include <getopt.h>
#include "dlib/bam_util.h"

int usage(char **argv, int retcode=EXIT_FAILURE) {
    fprintf(stderr, "%s <-l output_compression_level> in.bam out.bam\n"
                    "Use - for stdin or stdout.\n", argv[0]);
    return retcode;
}

int test_fail(bam1_t *b1, bam1_t *b2, void *data) {
    assert(bam_itag(b1, "FP") == 0 ? b1->core.flag & BAM_FQCFAIL: 1);
    assert(bam_itag(b2, "FP") == 0 ? b2->core.flag & BAM_FQCFAIL: 1);
    return 0;
}

int main(int argc, char *argv[]) {
    if(argc < 3) {
        return usage(argv);
    }
    if(strcmp(argv[1], "--help") == 0) {
        return usage(argv, EXIT_SUCCESS);
    }
    int c;
    char out_mode[4] = "wb";
    while((c = getopt(argc, argv, "l:h?")) > -1) {
        switch(c) {
        case 'l': out_mode[2] = atoi(optarg) % 10 + '0'; break;
        case 'h': case '?': return usage(argv, EXIT_SUCCESS);
        }
    }
    if(argc - 2 != optind) {
        LOG_EXIT("Required: precisely two positional arguments (in bam, out bam).\n");
    }
    // Actually this function. You can't really apply a null function....
    std::function<int (bam1_t *, void *)> fn = NULL;
    // Actually create your type for data and then provide it if needed.
    void *data = NULL;
    dlib::BamHandle inHandle(argv[optind]);
    dlib::BamHandle outHandle(argv[optind + 1], inHandle.header, "wb");
    dlib::abstract_pair_iter(inHandle.fp, inHandle.header, outHandle.fp,
                             &test_fail, nullptr);
    return EXIT_SUCCESS;
}
