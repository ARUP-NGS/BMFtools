/*  unclipped.c -- a postprocessor for bwa to facilitate MEI calls.
*/

#include "bam_util.h"

void add_unclipped(samFile *in, bam_hdr_t *hdr, samFile *ofp)
{
	abstract_pair_iter(in, hdr, ofp, &add_unclipped_mate_starts);
}

static int unclipped_usage(char *argv[]) {
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage: bmftools %s <input.namesrt.bam> <output.bam>\n\n", argv[0]);
	fprintf(stderr, "Opts:\n-l	 Sets bam compression level. (Valid: 1-9).\n");
	fprintf(stderr, "Set output.bam to \'-\' or \'stdout\' to pipe results.\n");
	fprintf(stderr, "Set input.amesrt.bam to \'-\' or \'stdoin\' to read from stdin.\n");

	sam_global_opt_help(stderr, "-....");
	exit(EXIT_FAILURE);
}

typedef struct option option_t;


int mark_unclipped_main(int argc, char *argv[])
{
	int c;
	samFile *in, *out;
	bam_hdr_t *header;
	char wmode[3] = {'w', 'b', '\0'};
	sam_global_args ga;
	memset(&ga, 0, sizeof(ga));

	static option_t lopts[] = {
		SAM_OPT_GLOBAL_OPTIONS('-', 0, 0, 0, 0),
		{ NULL, 0, NULL, 0 }
	};

	int level;
	while ((c = getopt_long(argc, argv, "l:?h", lopts, NULL)) >= 0) {
		switch (c) {
		case 'l': level = atoi(optarg); wmode[2] = level + '0'; break;
		default:  if (parse_sam_global_opt(c, optarg, lopts, &ga) == 0) break;
			/* else fall-through */
		case 'h': // Fall-through
		case '?': return unclipped_usage(argv);
		}
	}
	if (optind + 2 > argc)
		unclipped_usage(argv);

	in = sam_open_format(argv[optind], "r", &ga.in);
	header = sam_hdr_read(in);
	if (header == NULL || header->n_targets == 0) {
		fprintf(stderr, "[E:%s] input SAM does not have header. Abort!\n", __func__);
		return 1;
	}

	sam_open_mode(wmode+1, argv[optind+1], NULL);
	out = sam_open_format(argv[optind+1], wmode, &ga.out);
	if (in == 0 || out == 0) {
		fprintf(stderr, "[E:%s] fail to read/write input files\n", __func__);
		return 1;
	}
	sam_hdr_write(out, header);

	add_unclipped(in, header, out);
	bam_hdr_destroy(header);
	if(sam_close(in)) {
		fprintf(stderr, "[E:%s] Failed to close input file. Abort mission!\n", __func__);
		exit(EXIT_FAILURE);
	}
	if(sam_close(out)) {
		fprintf(stderr, "[E:%s] Failed to close output file. Abort mission!\n", __func__);
		exit(EXIT_FAILURE);
	}
	return 0;
}
