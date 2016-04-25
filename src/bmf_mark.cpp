#include <getopt.h>
#include "dlib/bam_util.h"

namespace BMF {

    namespace {
        struct mark_settings_t {
            // I might add more options later, hence the use of the bitfield.
            uint32_t remove_qcfail:1;
            uint32_t min_insert_length:8;
            double min_frac_unambiguous;
            mark_settings_t() :
                remove_qcfail(0),
                min_insert_length(0),
                min_frac_unambiguous(0.0)
            {
            }
        };

    }

    static inline int add_multiple_tags(bam1_t *b1, bam1_t *b2, void *data)
    {
        if(UNLIKELY(strcmp(bam_get_qname(b1), bam_get_qname(b2))))
            LOG_EXIT("Is this bam namesorted? These reads have different names.\n");
        int ret = 0;
        dlib::add_unclipped_mate_starts(b1, b2);
        dlib::add_mate_SA_tag(b1, b2);
        dlib::add_qseq_len(b1, b2);
        dlib::add_fraction_aligned(b1, b2);
        // Fails the reads if remove_qcfail is set and bitseq_qcfail returns 1
        ret |= (dlib::bitset_qcfail(b1, b2) && ((mark_settings_t *)data)->remove_qcfail);
        if(((mark_settings_t *)data)->min_insert_length)
            if(b1->core.isize)
                ret |= std::abs(b1->core.isize) < ((mark_settings_t *)data)->min_insert_length;
        ret |= dlib::filter_n_frac(b1, b2, ((mark_settings_t *)data)->min_frac_unambiguous);
#if !NDEBUG
        if(((mark_settings_t *)data)->remove_qcfail) {
            if(bam_itag(b1, "FP") == 0) {
                assert(ret);
            }
        }
#endif
        return ret;
    }

    static int mark_usage() {
        fprintf(stderr,
                        "Adds positional bam tags for a read and its mate for bmftools rsq and bmftools infer.\n"
                        "Meant primarily for piping to avoid I/O. Default compression is therefore 0. Typical compression for writing to disk: 6.\n"
                        "\tSU: Self Unclipped start.\n"
                        "\tMU: Mate Unclipped start.\n"
                        "\tms: Mate SA Tag. Supplemental Alignment tag. (Only if the mate has an SA tag).\n"
                        "Required for bmftools rsq using unclipped start.\n"
                        "Required for bmftools infer.\n"
                        "Usage: bmftools mark <opts> <input.namesrt.bam> <output.bam>\n\n"
                        "Flags:\n-l     Sets bam compression level. (Valid: 1-9). Default: 0.\n"
                        "-q    Skip read pairs which fail.\n"
                        "-d    Set bam compression level to default (6).\n"
                        "-i    Skip read pairs whose insert size is less than <INT>.\n"
                        "-u    Skip read pairs where both reads have a fraction of unambiguous base calls >= <FLOAT>\n"
                        "Set input.namesrt.bam to \'-\' or \'stdin\' to read from stdin.\n"
                        "Set output.bam to \'-\' or \'stdout\' or omit to stdout.\n"
                );
        exit(EXIT_FAILURE);
    }

    int mark_main(int argc, char *argv[])
    {
        char wmode[4]{"wb0"};
        int c;
        mark_settings_t settings;
        while ((c = getopt(argc, argv, "l:i:u:dq?h")) >= 0) {
            switch (c) {
            case 'u':
                settings.min_frac_unambiguous = atof(optarg); break;
            case 'q':
                settings.remove_qcfail = 1; break;
            case 'i':
                settings.min_insert_length = (uint32_t)atoi(optarg); break;
            case 'l':
                wmode[2] = atoi(optarg)%10 + '0';
                LOG_DEBUG("Now emitting output with compression level %c.\n", wmode[2]);
                break;
            case 'd':
                sprintf(wmode, "wb");
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

        dlib::BamHandle inHandle(in);
        dlib::BamHandle outHandle(out, inHandle.header, "wb");
        dlib::abstract_pair_iter(inHandle.fp, inHandle.header, outHandle.fp,
                                 &add_multiple_tags, &settings);

        LOG_INFO("Successfully complete bmftools mark.\n");
        return EXIT_SUCCESS;
    }

}
