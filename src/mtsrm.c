/*
 * Multi-threaded sort/rescue/merge
 *
 * Calls bmfsort on an input file or stdin, splits
 */

char **tmpfnames(char *split_prefix, bam_hdr_t *hdr) {
	char **ret = (char **)calloc(hdr->n_targets + 1, sizeof(char *));
	char buf[200];
	int i;
	for(i = 0; i < hdr->n_targets; ++i) {
		sprintf(buf, "%s.%s.bam", split_prefix, hdr->target_name[i]);
		ret[i] = strdup(buf);
	}
	sprintf(buf, "%s.unmapped.bam", split_prefix, hdr->target_name[i]);
	ret[i] = strdup(buf);
	return ret;
}

void srm_usage(FILE *fp, int retcode)
{
	fprintf(fp, "srm usage not written. Eh.\n");
	exit(retcode);
}


typedef struct {
	int mmlim;
	int is_se;
	int realign_unchanged;
	int cmpkey;
	int compression;
	char wmode[3];
	sam_global_args ga;
	bam_hdr_t *hdr;
	char **infnames;
	char **outfnames;
} pr_dispatcher_t;

int bam_pr_dispatcher(char *infname, char *outfname, int mmlim, int cmpkey, int is_se, int realign_unchanged,
					  int compression) {
    pr_dispatcher_t *dispatcher = (pr_dispatcher_t *)calloc(1, sizeof(pr_dispatcher_t));
    dispatcher->wmode = "wb";
    dispatcher->mmlim = mmlim;
    dispatcher->cmpkey = cmpkey;
    dispatcher->is_se = is_se;
    dispatcher->realign_unchanged = realign_unchanged;
    dispatcher->compression = compression;
    dispatcher->ga = SAM_GLOBAL_ARGS_INIT;
    samFile *tmpfp = sam_open_format(infname, "rb", &dispatcher->ga.in);
    dispatcher->hdr = sam_hdr_read(tmpfp);
}

int main(int argc, char *argv[]) {
	int c;
	char *inbam = NULL;
	char *final_bam = NULL;
	while ((c = getopt(argc, argv, "o:h?")) >= 0) {
		switch (c) {
		case 'o':
			final_bam = strdup(optarg);
			break;
		case '?':
		case 'h':
			srm_usage(stderr, EXIT_SUCCESS);
		default:
			srm_usage(stderr, EXIT_FAILURE);
		}
	}
	if(!final_bam) {
		fprintf(stderr, "Unset path for final bam. Writing to stdout.\n");
		final_bam = strdup("-");
	}
	cond_free(inbam);
	cond_free(final_bam);
	return 0;
}
