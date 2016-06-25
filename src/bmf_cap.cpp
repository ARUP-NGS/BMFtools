#include <getopt.h>
#include "dlib/bam_util.h"

namespace bmf {

static int cap_usage() {
    fprintf(stderr,
                    "\nCaps quality scores from PV tags to facilitate working with barcode-agnostic tools.\n"
                    "Values >= cap for PV will have the quality"
                    "strings set to the cap value. All others are set to 2 (#).\n"
                    "Usage: bmftools cap <input.bam> <output.bam>\n\n"
                    "Opts:\n-l     Sets bam compression level. (Valid: 1-9).\n"
                    "-c: set PV cap [uint32_t]. (Passing base qualities set to 93, < set to 2).\n"
                    "-m: set minFM to pass reads. [uint32_t].\n"
                    "-f: set minimum fraction agreed. [double].\n"
                    "-t: set maximum permitted phred score. [int, coerced to char].\n"
                    "-d: Flag to use existing quality scores instead of setting all below a threshold to 2.\n"
                    "Set output.bam to \'-\' or \'stdout\' to pipe results.\n"
                    "Set input.csrt.bam to \'-\' or \'stdin\' to read from stdin.\n"
            );
    return EXIT_FAILURE;
}

struct cap_settings_t {
    int32_t minFM;
    uint32_t minPV;
    double minFrac;
    int cap;
    int dnd;  // Do not disturb. Use uncapped unless >= minPV.
};


static inline int cap_bam_dnd(bam1_t *b, cap_settings_t *settings) {
    int i;
    uint32_t *PV((uint32_t *)dlib::array_tag(b, "PV"));
    uint32_t *FA((uint32_t *)dlib::array_tag(b, "FA"));
    const int FM(bam_itag(b, "FM"));
    if(FM < settings->minFM)
        return 1;
    const int l_qseq (b->core.l_qseq);
    char *qual((char *)bam_get_qual(b));
    if(b->core.flag & BAM_FREVERSE) {
        for(i = 0; i < l_qseq; ++i) {
            if(PV[i] >= settings->minPV && static_cast<double>(FA[i]) / FM >= settings->minFrac) {
                qual[l_qseq - 1 - i] = settings->cap;
            }
        }
    } else {
        for(i = 0; i < l_qseq; ++i)
            if(PV[i] >= settings->minPV && static_cast<double>(FA[i]) / FM >= settings->minFrac)
                qual[i] = settings->cap;
    }
    return 0;
}


static inline int cap_bam_q(bam1_t *b, cap_settings_t *settings) {
    uint32_t *const PV((uint32_t *)dlib::array_tag(b, "PV"));
    uint32_t *const FA((uint32_t *)dlib::array_tag(b, "FA"));
    const int FM(bam_itag(b, "FM"));
    if(FM < (int)settings->minFM)
        return 1;
    char *qual((char *)bam_get_qual(b));
    int i(0);
    // If the thresholds are failed, set the quality to 2 to make them
    if(b->core.flag & BAM_FREVERSE) {
        int l_qseq(b->core.l_qseq);
        for(; l_qseq >= 0; ++i) {
            if(PV[i] >= settings->minPV && static_cast<double>(FA[i]) / FM >= settings->minFrac)
                qual[--l_qseq] = settings->cap;
            else
                qual[--l_qseq] = 2;
        }
    } else {
        for(; i < b->core.l_qseq; ++i) {
            if(PV[i] >= settings->minPV && static_cast<double>(FA[i]) / FM >= settings->minFrac) 
                qual[i] = settings->cap;
            else
                qual[i] = 2;
        }
    }
    return 0;
}


int cap_main(int argc, char *argv[])
{
    cap_settings_t settings{0};
    settings.cap = 93;
    int c;
    char wmode[4]{"wb"};

    int level(6);
    while ((c = getopt(argc, argv, "t:f:m:l:c:dh?")) >= 0) {
        switch (c) {
        // mod 10 to handle a user error of negative or excessively high compression level.
        case 'l':
            level = atoi(optarg) % 10; wmode[2] = level + '0';
            LOG_DEBUG("Emitting bam with compression level %i.\n", level);
            break;
        case 'm': settings.minFM = strtoul(optarg, nullptr, 10); break;
        case 'c': settings.minPV = strtoul(optarg, nullptr, 10); break;
        case 'f': settings.minFrac = atof(optarg); break;
        case 't':
            settings.cap = (char)atoi(optarg);
            if(settings.cap > 93) {
                LOG_EXIT("Requested cap %i > maximum cap score 93. Abort!\n", settings.cap);
            }
            break;
        case 'd': settings.dnd = 1; break;
        case 'h': case '?': cap_usage(); return EXIT_SUCCESS;
        }
    }
    if (optind + 2 > argc)
        return cap_usage();

    if(settings.minPV == 0 && settings.minFM == 0 && settings.minFrac == 0.0) {
        fprintf(stderr, "[E:%s] All caps cannot be set to 0 (default values). [Required parameter] See usage.\n", __func__);
        return cap_usage();
    }
    int ret(dlib::bam_apply_function(argv[optind], argv[optind+1],
                                     settings.dnd ? (single_aux_fn)&cap_bam_dnd
                                                  : (single_aux_fn)&cap_bam_q,
                                     (void *)&settings, wmode));
    if(ret) {
    	fprintf(stderr, "bmftools cap returned non-zero exit status %i.\n", ret);
    	return ret;
    }
    LOG_INFO("Successfully completed bmftools cap!\n");
    return ret;
}

} /* namespace bmf */
