/*  unclipped.c -- a postprocessor for bwa to facilitate MEI calls.
*/

#include "dlib/bam_util.h"
#include "include/sam_opts.h"

static inline int add_multiple_tags(bam1_t *b1, bam1_t *b2, void *data)
{
    add_unclipped_mate_starts(b1, b2);
    add_sc_lens(b1, b2);
    add_fraction_aligned(b1, b2);
    return 0;
}

static int unclipped_usage(char *argv[]) {
    fprintf(stderr, "\n");
    fprintf(stderr, "Usage: bmftools %s <input.namesrt.bam> <output.bam>\n\n", argv[0]);
    fprintf(stderr, "Opts:\n-l     Sets bam compression level. (Valid: 1-9).\n");
    fprintf(stderr, "Set output.bam to \'-\' or \'stdout\' to pipe results.\n");
    fprintf(stderr, "Set input.namesrt.bam to \'-\' or \'stdin\' to read from stdin.\n");

    exit(EXIT_FAILURE);
}

int mark_unclipped_main(int argc, char *argv[])
{
    char wmode[4] = {'w', 'b', '\0', '\0'};
    int c;
    while ((c = getopt(argc, argv, "l:?h")) >= 0) {
        switch (c) {
        case 'l': wmode[2] = atoi(optarg)%10 + '0'; break;
        case 'h': case '?': unclipped_usage(argv);
        }
    }

    if (optind + 2 > argc) unclipped_usage(argv);

    return dlib::bam_pair_apply_function(argv[optind], argv[optind+1],
                                         add_multiple_tags, NULL, wmode);
    return EXIT_SUCCESS;
}
