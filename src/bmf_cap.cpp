/*  cap_qscore.cpp -- a postprocessor for bmftools dmp which
 *  prepares our sophisticated error model for downstream analysis
 *  for BMF-agnostic tools.
*/
#include "bmf_cap.h"



	void cap_qscore_core(samFile *in, bam_hdr_t *hdr, samFile *out, cap_settings *settings)
	{
		abstract_single_filter(in, hdr, out, settings->dnd ? (single_aux_check)&cap_bam_dnd: (single_aux_check)&cap_bam_q, (void *)settings);
	}


	static int cap_qscore_usage(char *argv[]) {
		fprintf(stderr, "\ncap_qscore takes a bam and caps quality scores from PV tags "
				"to facilitate working with outside tools.\nValues >= cap for PV will have the quality"
				"strings set to the cap value. All others are set to 2 (#).\n");
		fprintf(stderr, "Usage:  %s <input.bam> <output.bam>\n\n", argv[0]);
		fprintf(stderr, "Opts:\n-l	 Sets bam compression level. (Valid: 1-9).\n");
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
	samFile *in, *out;
	bam_hdr_t *header;
	char wmode[3] = {'w', 'b', 0};
	sam_global_args ga;
	memset(&ga, 0, sizeof(ga));

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
	if(!(in = sam_open_format(argv[optind], "r", &ga.in))) {
		fprintf(stderr, "[E:%s] input bam '%s' could not be read. Abort!\n", __func__, argv[optind]);
		return EXIT_FAILURE;
	}
	if ((header = sam_hdr_read(in)) == NULL || header->n_targets == 0) {
		fprintf(stderr, "[E:%s] input SAM does not have header. Abort!\n", __func__);
		return EXIT_FAILURE;
	}

	sam_open_mode(wmode+1, argv[optind+1], NULL);
	out = sam_open_format(argv[optind+1], wmode, &ga.out);
	if (in == 0 || out == 0) {
		fprintf(stderr, "[pair_mark] fail to read/write input files\n");
		return EXIT_FAILURE;
	}
	sam_hdr_write(out, header);

	cap_qscore_core(in, header, out, &settings);
	bam_hdr_destroy(header);
	sam_close(in); sam_close(out);
	return 0;
}
