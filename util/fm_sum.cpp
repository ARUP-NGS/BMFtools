#include "dlib/bam_util.h"
#include <getopt.h>

int usage(char **argv, int retcode=EXIT_FAILURE) {
    fprintf(stderr, "%s <-l output_compression_level> in.bam out.bam\n"
                    "Use - for stdin or stdout.\n", argv[0]);
    return retcode;
}

int main(int argc, char *argv[]) {
    if(argc < 2)
        return usage(argv);
    if(strcmp(argv[1], "--help") == 0) {
        return usage(argv, EXIT_SUCCESS);
    }
    int c;
    while((c = getopt(argc, argv, "l:h?")) > -1) {
        switch(c) {
        case 'h': case '?': return usage(argv, EXIT_SUCCESS);
        }
    }
    if(argc - 2 != optind)
        LOG_EXIT("Required: precisely two positional arguments (in bam, out bam).\n");
    dlib::BamHandle in(argv[1]);
    bam1_t *b(bam_init1());
    size_t count(0);
    while(sam_read1(in.fp, in.header, b) >= 0) count += bam_itag(b,"FM");
    fprintf(stdout, "%s: %lu\n", argv[1], count);
    bam_destroy1(b);
    return EXIT_SUCCESS;
}
