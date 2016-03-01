/*  cap_qscore.cpp -- a postprocessor for bmftools dmp which
 *  prepares our sophisticated error model for downstream analysis
 *  for BMF-agnostic tools.
*/
#include "bmf_cap.h"

static int cap_qscore_usage(char *argv[]) {
    fprintf(stderr, "\ncap_qscore takes a bam and caps quality scores from PV tags "
            "to facilitate working with outside tools.\nValues >= cap for PV will have the quality"
            "strings set to the cap value. All others are set to 2 (#).\n");
    fprintf(stderr, "Usage:  %s <input.bam> <output.bam>\n\n", argv[0]);
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


int cap_qscore_main(int argc, char *argv[])
{
    cap_settings settings = {0};
    settings.cap = 93;
    int c;
    char wmode[3] = {'w', 'b', 0};

    int level;
    while ((c = getopt(argc, argv, "t:f:m:l:c:dh?")) >= 0) {
        switch (c) {
        case 'l': level = atoi(optarg) % 10; wmode[2] = level + '0'; break;
        case 'm': settings.minFM = strtoul(optarg, NULL, 10); break;
        case 'c': settings.minPV = strtoul(optarg, NULL, 10); break;
        case 'f': settings.minFrac = atof(optarg); break;
        case 't':
            if(atoi(optarg) > 93) {
                fprintf(stderr, "Hey, this qscore is too high. 93 max!\n");
                exit(EXIT_FAILURE);
            }
            settings.cap = (char)atoi(optarg); break;
        case 'd': settings.dnd = 1; break;
        case 'h': // fall-through
        case '?': return cap_qscore_usage(argv);
        }
    }
    if (optind + 2 > argc)
        return cap_qscore_usage(argv);

    if(settings.minPV == 0 && settings.minFM == 0 && settings.minFrac == 0.0) {
        fprintf(stderr, "[E:%s] All caps cannot be set to 0 (default value). [Required parameter] See usage.\n", __func__);
        return cap_qscore_usage(argv);
    }
    int ret = dlib::bam_apply_function(argv[optind], argv[optind+1],
                                       settings.dnd ? (single_aux_check)&cap_bam_dnd
                                                    : (single_aux_check)&cap_bam_q,
                                       (void *)&settings, wmode);
    return ret;
}
