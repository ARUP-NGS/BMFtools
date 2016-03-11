#include "dlib/bam_util.h"

namespace BMF {

    namespace {
        struct mark_settings_t {
            // I might add more options later, hence the use of the bitfield.
            uint32_t remove_qcfail:1;
            uint32_t min_insert_length:8;
        };
    }

    static inline int add_multiple_tags(bam1_t *b1, bam1_t *b2, void *data)
    {
        int ret = 0;
        dlib::add_unclipped_mate_starts(b1, b2);
        dlib::add_sc_lens(b1, b2);
        dlib::add_fraction_aligned(b1, b2);
        dlib::add_qseq_len(b1, b2);
/*
        if(((mark_settings_t *)data)->remove_qcfail) {
            ret |= dlib::bitset_qcfail(b1, b2);
        } else {
            dlib::bitset_qcfail(b1, b2);
        }
*/
        if(((mark_settings_t *)data)->min_insert_length) {
            if(b1->core.isize && // Non-zero insert size
               std::abs(b1->core.isize) < ((mark_settings_t *)data)->min_insert_length)
                ret = 1; // Fail it!
        }
        return ret;
    }

    static int mark_usage() {
        fprintf(stderr,
                        "Adds positional bam tags for a read and its mate for bmftools rsq and bmftools infer.\n"
                        "\tSU: Self Unclipped start.\n"
                        "\tMU: Mate Unclipped start.\n"
                        "\tLM: Mate length.\n"
                        "Required for bmftools rsq using unclipped start.\n"
                        "Required for bmftools infer.\n"
                        "Usage: bmftools mark <opts> <input.namesrt.bam> <output.bam>\n\n"
                        "Flags:\n-l     Sets bam compression level. (Valid: 1-9).\n"
                        "-q    Skip read pairs which fail.\n"
                        "-i    Skip read pairs whose insert size is less than this number.\n"
                        "Set input.namesrt.bam to \'-\' or \'stdin\' to read from stdin.\n"
                        "Set output.bam to \'-\' or \'stdout\' or omit to stdout.\n"
                );
        exit(EXIT_FAILURE);
    }

    int mark_main(int argc, char *argv[])
    {
        char wmode[4]{"wb"};
        int c;
        mark_settings_t settings{0};
        while ((c = getopt(argc, argv, "l:i:q?h")) >= 0) {
            switch (c) {
            case 'q':
                settings.remove_qcfail = 1; break;
            case 'i':
                settings.min_insert_length = (uint32_t)atoi(optarg); break;
            case 'l':
                    wmode[2] = atoi(optarg)%10 + '0';
                    LOG_DEBUG("Now emitting output with compression level %c.\n", wmode[2]);
                    break;
            case '?': case 'h': mark_usage(); // Exits. No need for a break.
            }
        }

        char *in = (char *)"-", *out = (char *)"-";
        if(optind + 2 == argc) {
            in = argv[optind];
            out = argv[optind + 1];
        } else if(optind + 1 == argc) {
            LOG_INFO("No outfile provided. Defaulting to stdout.\n");
            in = argv[optind];
        } else {
            LOG_INFO("No input or output bam provided! Defaulting stdin and stdout.\n");
        }

        int ret = dlib::bam_pair_apply_function(in, out,
                                                add_multiple_tags, (void *)&settings, wmode);
        LOG_INFO("Successfully complete bmftools mark.\n");
        return ret;
    }

}
