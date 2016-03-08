#include "dlib/bam_util.h"
#include "include/sam_opts.h"

namespace BMF {

static inline int add_multiple_tags(bam1_t *b1, bam1_t *b2, void *data)
{
    dlib::add_unclipped_mate_starts(b1, b2);
    dlib::add_sc_lens(b1, b2);
    dlib::add_fraction_aligned(b1, b2);
    return 0;
}

static int unclipped_usage() {
    fprintf(stderr,
                    "Adds the unclipped start position for each read and its mate as tags:\n"
                    "\tSU: Self Unclipped start.\n"
                    "\tMU: Mate Unclipped start.\n"
                    "Required for bmftools rsq using unclipped start.\n"
                    "Usage: bmftools mark_unclipped <opts> <input.namesrt.bam> <output.bam>\n\n"
                    "Flags:\n-l     Sets bam compression level. (Valid: 1-9).\n"
                    "Set output.bam to \'-\' or \'stdout\' to pipe results.\n"
                    "Set input.namesrt.bam to \'-\' or \'stdin\' to read from stdin.\n"
            );
    exit(EXIT_FAILURE);
}

int mark_unclipped_main(int argc, char *argv[])
{
    char wmode[4] = "wb";
    int c;
    while ((c = getopt(argc, argv, "l:?h")) >= 0) {
        switch (c) {
        case 'l':
                wmode[2] = atoi(optarg)%10 + '0';
                LOG_DEBUG("Now emitting output with compressio level %c.\n", wmode[2]);
                break;
        case '?': case 'h': unclipped_usage(); break;
        }
    }

    if (optind + 2 > argc) unclipped_usage();

    int ret = dlib::bam_pair_apply_function(argv[optind], argv[optind+1],
                                            add_multiple_tags, NULL, wmode);
    LOG_INFO("Successfully complete bmftools mark_unclipped.\n");
    return ret;
}

}
