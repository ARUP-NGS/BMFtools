/*  cap_qscore.cpp -- a postprocessor for bmftools dmp which
 *  prepares our sophisticated error model for downstream analysis
 *  for BMF-agnostic tools.
*/
#include "dlib/bam_util.h"
#include "include/sam_opts.h"

namespace BMF {

    static int cap_qscore_usage() {
        fprintf(stderr, "\nTakes a bam and caps quality scores from PV tags "
                "to facilitate working with outside tools.\nValues >= cap for PV will have the quality"
                "strings set to the cap value. All others are set to 2 (#).\n");
        fprintf(stderr, "Usage: bmftools cap <input.bam> <output.bam>\n\n");
        fprintf(stderr, "Opts:\n-l     Sets bam compression level. (Valid: 1-9).\n");
        fprintf(stderr, "-c: set PV cap [uint32_t]. (Passing base qualities set to 93, < set to 2).\n");
        fprintf(stderr, "-m: set minFM to pass reads. [uint32_t].\n");
        fprintf(stderr, "-f: set minimum fraction agreed. [double].\n");
        fprintf(stderr, "-t: set maximum permitted phred score. [int, coerced to char].\n");
        fprintf(stderr, "-d: Flag to use existing quality scores instead of setting all below a threshold to 2.\n");
        fprintf(stderr, "Set output.bam to \'-\' or \'stdout\' to pipe results.\n");
        fprintf(stderr, "Set input.csrt.bam to \'-\' or \'stdin\' to read from stdin.\n");
        return EXIT_FAILURE;
    }

    struct cap_settings_t {
        int32_t minFM;
        uint32_t minPV;
        double minFrac;
        int cap;
        int dnd; // Do not disturb. Use uncapped unless greater than given minPV.
    };


    static inline int cap_bam_dnd(bam1_t *b, cap_settings_t *settings) {
        int i;
        uint32_t *PV = (uint32_t *)dlib::array_tag(b, "PV");
        uint32_t *FA = (uint32_t *)dlib::array_tag(b, "FA");
        const int FM = bam_itag(b, "FM");
        if(FM < settings->minFM)
            return 1;
        const int l_qseq = b->core.l_qseq;
        char *qual = (char *)bam_get_qual(b);
        if(b->core.flag & BAM_FREVERSE)
            for(i = 0; i < l_qseq; ++i)
                if(PV[i] >= settings->minPV && (double)FA[i] / FM >= settings->minFrac)
                    qual[l_qseq - 1 - i] = settings->cap;
        else
            for(i = 0; i < l_qseq; ++i)
                if(PV[i] >= settings->minPV && (double)FA[i] / FM >= settings->minFrac)
                    qual[i] = settings->cap;
        return 0;
    }

    static inline int cap_bam_q(bam1_t *b, cap_settings_t *settings) {
        uint32_t *const PV = (uint32_t *)dlib::array_tag(b, "PV");
        uint32_t *const FA = (uint32_t *)dlib::array_tag(b, "FA");
        const int FM = bam_itag(b, "FM");
        if(FM < (int)settings->minFM)
            return 1;
        char *qual = (char *)bam_get_qual(b);
        int i;
        const int l_qseq = b->core.l_qseq;
        // If the thresholds are failed, set the quality to 2 to make them
        if(b->core.flag & BAM_FREVERSE)
            for(i = 0; i < l_qseq; ++i)
                qual[l_qseq - 1 - i] = (PV[i] >= settings->minPV && (double)FA[i] / FM >= settings->minFrac) ? settings->cap: 2;
        else
            for(i = 0; i < l_qseq; ++i)
                qual[i] = (PV[i] >= settings->minPV && (double)FA[i] / FM >= settings->minFrac) ? settings->cap: 2;
        return 0;
    }



    int cap_qscore_main(int argc, char *argv[])
    {
        cap_settings_t settings = {0};
        settings.cap = 93;
        int c;
        char wmode[4] = "wb";

        int level;
        while ((c = getopt(argc, argv, "t:f:m:l:c:dh?")) >= 0) {
            switch (c) {
            // mod 10 to handle a user error of negative or excessively high compression level.
            case 'l':
                level = atoi(optarg) % 10; wmode[2] = level + '0';
                LOG_DEBUG("Emitting bam with compression level %i.\n", level);
                break;
            case 'm': settings.minFM = strtoul(optarg, NULL, 10); break;
            case 'c': settings.minPV = strtoul(optarg, NULL, 10); break;
            case 'f': settings.minFrac = atof(optarg); break;
            case 't':
                settings.cap = (char)atoi(optarg);
                if(settings.cap > 93) {
                    LOG_EXIT("Requested cap %i > maximum cap score 93. Abort!\n", settings.cap);
                }
                break;
            case 'd': settings.dnd = 1; break;
            case 'h': case '?': cap_qscore_usage(); return EXIT_SUCCESS;
            }
        }
        if (optind + 2 > argc)
            return cap_qscore_usage();

        if(settings.minPV == 0 && settings.minFM == 0 && settings.minFrac == 0.0) {
            fprintf(stderr, "[E:%s] All caps cannot be set to 0 (default value). [Required parameter] See usage.\n", __func__);
            return cap_qscore_usage();
        }
        return dlib::bam_apply_function(argv[optind], argv[optind+1],
                                        settings.dnd ? (single_aux_check)&cap_bam_dnd
                                                     : (single_aux_check)&cap_bam_q,
                                        (void *)&settings, wmode);
    }

}
