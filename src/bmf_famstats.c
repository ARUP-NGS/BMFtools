#include "bmf_famstats.h"

int RVWarn = 1;

int frac_usage_exit(FILE *fp, int code) {
	fprintf(fp, "bmftools famstats frac <opts> <in.bam>\n"
			"Opts:\n-m minFM to accept. REQUIRED.\n"
			"-h, -?: Return usage.\n");
	exit(code);
	return EXIT_FAILURE; // This never happens
}


static void print_hashstats(famstats_t *stats, FILE *fp)
{
	fprintf(fp, "#Family size\tNumber of families\n");
	for(stats->ki = kh_begin(stats->fm); stats->ki != kh_end(stats->fm); ++stats->ki) {
		if(!kh_exist(stats->fm, stats->ki)) continue;
		fprintf(fp, "%"PRIu64"\t%"PRIu64"\n", kh_key(stats->fm, stats->ki), kh_val(stats->fm, stats->ki));
	}
	fprintf(fp, "#RV'd in family\tNumber of families\n");
	for(stats->ki = kh_begin(stats->rc); stats->ki != kh_end(stats->rc); ++stats->ki) {
		if(!kh_exist(stats->rc, stats->ki)) continue;
		fprintf(fp, "%"PRIu64"\t%"PRIu64"\n", kh_key(stats->rc, stats->ki), kh_val(stats->rc, stats->ki));
	}
}


void target_usage_exit(FILE *fp, int success)
{
	fprintf(fp, "Usage: bmftools famstats target <opts> <in.bam>\nOpts:\n-b Path to bed file.\n"
			"-p padding. Number of bases around bed regions to pad. Default: 25.\n"
			"-h, -?: Return usage.\n");
	exit(success);
}

static inline void target_loop(bam1_t *b, khash_t(bed) *bed, uint64_t *fm_target, uint64_t *total_fm)
{
	if(!bam_aux2i(bam_aux_get(b, "FP"))) return;
	const int FM = bam_aux2i(bam_aux_get(b, "FM"));
	*total_fm += FM;
	if(bed_test(b, bed))
		*fm_target += FM;
}


int target_main(int argc, char *argv[])
{
	samFile *fp;
	bam_hdr_t *header;
	int c;
	khash_t(bed) *bed = NULL;
	char *bedpath = NULL;
	uint32_t padding = (uint32_t)-1;
	uint64_t notification_interval = 1000000;

	if(argc == 1) target_usage_exit(stderr, EXIT_SUCCESS);
	else if(argc < 4) {
		target_usage_exit(stderr, EXIT_FAILURE);
	}

	if(strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0) target_usage_exit(stderr, EXIT_SUCCESS);

	while ((c = getopt(argc, argv, "b:p:n:h?")) >= 0) {
		switch (c) {
		case 'b':
			bedpath = strdup(optarg);
			break;
		case 'p':
			padding = strtoul(optarg, NULL, 0);
			break;
		case 'n':
			notification_interval = strtoull(optarg, NULL, 0);
			break;
		case '?':
		case 'h':
			target_usage_exit(stderr, EXIT_SUCCESS);
		default:
			target_usage_exit(stderr, EXIT_FAILURE);
		}
	}


	if(padding == (uint32_t)-1) {
		padding = 25;
		fprintf(stderr, "[famstat_target_main]: Padding not set. Set to 25.\n");
	}

	if (argc != optind+1)
		(argc == optind) ? target_usage_exit(stdout, EXIT_SUCCESS): target_usage_exit(stderr, EXIT_FAILURE);

	if(!bedpath) {
		fprintf(stderr, "Bed path required for famstats target. See usage!\n");
		target_usage_exit(stderr, EXIT_FAILURE);
	}

	fp = sam_open(argv[optind], "r");
	if (fp == NULL) {
		fprintf(stderr, "[famstat_target_main]: Cannot open input file \"%s\"", argv[optind]);
		exit(EXIT_FAILURE);
	}

	header = sam_hdr_read(fp);
	if (header == NULL) {
		fprintf(stderr, "[famstat_target_main]: Failed to read header for \"%s\"\n", argv[optind]);
		exit(EXIT_FAILURE);
	}
	bed = parse_bed_hash(bedpath, header, padding);
	uint64_t fm_target = 0, total_fm = 0, count = 0;
	bam1_t *b = bam_init1();
	while (LIKELY((c = sam_read1(fp, header, b)) >= 0)) {
		target_loop(b, bed, &fm_target, &total_fm);
		if(UNLIKELY(++count % notification_interval == 0))
			fprintf(stderr, "[%s] Number of records processed: %lu.\n", __func__, count);
	}
	bam_destroy1(b);
	bam_hdr_destroy(header);
	sam_close(fp);
	bed_destroy_hash(bed);
	free(bedpath);
	fprintf(stdout, "Fraction of raw reads on target: %f. \nTotal raw reads: %"PRIu64". Raw reads on target: %"PRIu64".\n",
			(double)fm_target / total_fm, total_fm, fm_target);
	return 0;
}


static void print_stats(famstats_t *stats, FILE *fp)
{
	fprintf(fp, "#Number passing filters: %"PRIu64".\n", stats->n_pass);
	fprintf(fp, "#Number failing filters: %"PRIu64".\n", stats->n_fail);
	fprintf(fp, "#Summed FM (total founding reads): %"PRIu64".\n", stats->allfm_sum);
	fprintf(fp, "#Summed FM (total founding reads), (FM > 1): %"PRIu64".\n", stats->realfm_sum);
	fprintf(fp, "#Summed RV (total reverse-complemented reads): %"PRIu64".\n", stats->allrc_sum);
	fprintf(fp, "#Summed RV (total reverse-complemented reads), (FM > 1): %"PRIu64".\n", stats->realrc_sum);
	fprintf(fp, "#RV fraction for all read families: %lf.\n", (double)stats->allrc_sum / (double)stats->allfm_sum);
	fprintf(fp, "#RV fraction for real read families: %lf.\n", (double)stats->realrc_sum / (double)stats->realfm_sum);
	fprintf(fp, "#Mean Family Size (all)\t%lf\n", (double)stats->allfm_sum / (double)stats->allfm_counts);
	fprintf(fp, "#Mean Family Size (real)\t%lf\n", (double)stats->realfm_sum / (double)stats->realfm_counts);
	print_hashstats(stats, fp);
}

static inline void tag_test(const uint8_t *data, const char *tag)
{
	if(UNLIKELY(!data))
		fprintf(stderr, "Required bam tag '%s' not found. Abort mission!\n", tag),
		exit(EXIT_FAILURE);
}


static inline void famstat_loop(famstats_t *s, bam1_t *b, famstat_settings_t *settings)
{
	const uint8_t *data = bam_aux_get(b, "FM");
	//tag_test(data, "FM");
	const int FM = bam_aux2i(data);
	if(b->core.flag & 2944 || b->core.qual < settings->minMQ || FM < settings->minFM ||
		!bam_aux2i(bam_aux_get(b, "FP"))) {
		++s->n_fail;
		return;
		// 2944 is equivalent to BAM_FSECONDARY | BAM_FSUPPLEMENTARY | BAM_QCFAIL | BAM_FISREAD2
	}
	++s->n_pass;
	const uint8_t *rvdata = bam_aux_get(b, "RV");
	if(!rvdata && RVWarn) {
		RVWarn = 0;
		fprintf(stderr, "[%s]: Warning: RV tag not found. Continue.\n", __func__);
	}
	int RV = rvdata ? bam_aux2i(data): 0;

	if(FM > 1)
		++s->realfm_counts, s->realfm_sum += FM, s->realrc_sum += RV;
	++s->allfm_counts, s->allfm_sum += FM, s->allrc_sum += RV;

	if((s->ki = kh_get(fm, s->fm, FM)) == kh_end(s->fm))
		s->ki = kh_put(fm, s->fm, FM, &s->khr), kh_val(s->fm, s->ki) = 1;
	else
		++kh_val(s->fm, s->ki);
	if((s->ki = kh_get(rc, s->rc, RV)) == kh_end(s->rc))
		s->ki = kh_put(rc, s->rc, RV, &s->khr), kh_val(s->rc, s->ki) = 1;
	else
		++kh_val(s->rc, s->ki);
}


famstats_t *famstat_core(samFile *fp, bam_hdr_t *h, famstat_settings_t *settings)
{
	uint64_t count = 0;
	famstats_t *s;
	bam1_t *b;
	int ret;
	s = (famstats_t*)calloc(1, sizeof(famstats_t));
	s->fm = kh_init(fm);
	s->rc = kh_init(rc);
	s->data = NULL;
	b = bam_init1();
	while (LIKELY((ret = sam_read1(fp, h, b)) >= 0)) {
		famstat_loop(s, b, settings);
		if(UNLIKELY(++count % settings->notification_interval == 0))
			fprintf(stderr, "[%s] Number of records processed: %lu.\n", __func__, count);
	}
	bam_destroy1(b);
	if (ret != -1)
		fprintf(stderr, "[%s] Truncated file? Continue anyway.\n", __func__);
	return s;
}

static void usage_exit(FILE *fp, int exit_status)
{
	fprintf(fp, "Usage: bmftools famstats\n");
	fprintf(fp, "Subcommands: \nfm\tFamily Size stats\n"
			"frac\tFraction of raw reads in family sizes >= minFM parameter.\n"
			"target\tFraction of raw reads on target.\n");
	exit(exit_status);
}

static void fm_usage_exit(FILE *fp, int exit_status)
{
	fprintf(fp, "Usage: bmftools famstats fm <opts> <in.bam>\n");
	fprintf(fp, "-m Set minimum mapping quality. Default: 0.\n");
	fprintf(fp, "-f Set minimum family size. Default: 0.\n");
	exit(exit_status);
}

int fm_main(int argc, char *argv[])
{
	samFile *fp;
	bam_hdr_t *header;
	famstats_t *s;
	int c;


	famstat_settings_t *settings = (famstat_settings_t *)calloc(1, sizeof(famstat_settings_t));
	settings->minMQ = 0;
	settings->minFM = 0;
	settings->notification_interval = 1000000uL;

	while ((c = getopt(argc, argv, "m:f:n:h")) >= 0) {
		switch (c) {
		case 'm':
			settings->minMQ = atoi(optarg); break;
			break;
		case 'f':
			settings->minFM = atoi(optarg); break;
			break;
		case 'n':
			settings->notification_interval = strtoull(optarg, NULL, 0);
		case '?': // Fall-through!
		case 'h':
			fm_usage_exit(stderr, EXIT_SUCCESS);
		default:
			fm_usage_exit(stderr, EXIT_FAILURE);
		}
	}

	if (argc != optind+1) {
		if (argc == optind) fm_usage_exit(stdout, EXIT_SUCCESS);
		else fm_usage_exit(stderr, EXIT_FAILURE);
	}

	fprintf(stderr, "[famstat_main]: Running main with minMQ %i and minFM %i.\n", settings->minMQ, settings->minFM);

	fp = sam_open(argv[optind], "r");
	if (fp == NULL) {
		fprintf(stderr, "[famstat_main]: Cannot open input file \"%s\"", argv[optind]);
		exit(EXIT_FAILURE);
	}

	header = sam_hdr_read(fp);
	if (header == NULL) {
		fprintf(stderr, "[famstat_main]: Failed to read header for \"%s\"\n", argv[optind]);
		exit(EXIT_FAILURE);
	}
	s = famstat_core(fp, header, settings);
	print_stats(s, stdout);
	kh_destroy(fm, s->fm);
	kh_destroy(rc, s->rc);
	free(s);
	free(settings);
	bam_hdr_destroy(header);
	sam_close(fp);
	return 0;
}

static inline void frac_loop(bam1_t *b, int minFM, uint64_t *fm_above, uint64_t *fm_total)
{
	if((b->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY | BAM_FQCFAIL | BAM_FREAD2)) || bam_aux2i(bam_aux_get(b, "FP")) == 0)
		return;
	const uint8_t *const data = bam_aux_get(b, "FM");
	//tag_test(data, "FM");
	const int FM = bam_aux2i(data);
	*fm_total += FM;
	if(FM >= minFM) *fm_above += FM;
}


int frac_main(int argc, char *argv[])
{
	samFile *fp;
	bam_hdr_t *header;
	int c;
	uint32_t minFM = 0;
	uint64_t notification_interval = 1000000;

	if(argc < 4) frac_usage_exit(stderr, EXIT_FAILURE);
	if(strcmp(argv[1], "--help") == 0) frac_usage_exit(stderr, EXIT_SUCCESS);

	while ((c = getopt(argc, argv, "n:m:h?")) >= 0) {
		switch (c) {
		case 'm':
			minFM = (uint32_t)atoi(optarg); break;
			break;
		case 'n':
			notification_interval = strtoull(optarg, NULL, 0); break;
		case '?':
		case 'h':
			return frac_usage_exit(stderr, EXIT_SUCCESS);
		default:
			return frac_usage_exit(stderr, EXIT_FAILURE);
		}
	}

	if(!minFM) {
		fprintf(stderr, "minFM not set. frac_main meaningless without it. Result: 1.0.\n");
		return EXIT_FAILURE;
	}
	fprintf(stderr, "[famstat_frac_main]: Running frac main minFM %i.\n", minFM);

	if (argc != optind+1) {
		if (argc == optind) frac_usage_exit(stdout, EXIT_SUCCESS);
		else frac_usage_exit(stderr, EXIT_FAILURE);
	}
	fp = sam_open(argv[optind], "r");
	if (fp == NULL) {
		fprintf(stderr, "[famstat_frac_main]: Cannot open input file \"%s\"", argv[optind]);
		exit(EXIT_FAILURE);
	}

	header = sam_hdr_read(fp);
	if (header == NULL) {
		fprintf(stderr, "[famstat_main]: Failed to read header for \"%s\"\n", argv[optind]);
		exit(EXIT_FAILURE);
	}
	uint64_t fm_above = 0, total_fm = 0, count = 0;
	bam1_t *b = bam_init1();
	// Check to see if the required tags are present before starting.
	if((c = sam_read1(fp, header, b)) < 0) {
		fprintf(stderr, "[E:%s] Could not read initial record from input file '%s'. Abort!\n", __func__, argv[optind]);
		exit(EXIT_FAILURE);
	}
	check_bam_tag(b, "FP");
	check_bam_tag(b, "RV");
	check_bam_tag(b, "FM");
	check_bam_tag(b, "FA");
	frac_loop(b, minFM, &fm_above, &total_fm), ++count;
	while (LIKELY(c = sam_read1(fp, header, b)) >= 0) {
		frac_loop(b, minFM, &fm_above, &total_fm);
		if(UNLIKELY(!(++count % notification_interval)))
			fprintf(stderr, "[famstat_frac_core] Number of records processed: %lu.\n", count);
	}
	fprintf(stderr, "#Fraction of raw reads with >= minFM %i: %f.\n", minFM, (double)fm_above / total_fm);
	bam_destroy1(b);
	bam_hdr_destroy(header);
	sam_close(fp);
	return 0;
}

int famstats_main(int argc, char *argv[])
{
	if(argc < 2) {
		usage_exit(stderr, 1);
	}
	if(strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0) {
		usage_exit(stderr, 0);
	}
	if(strcmp(argv[1], "fm") == 0) {
		return fm_main(argc - 1, argv + 1);
	}
	if(strcmp(argv[1], "frac") == 0) {
		return frac_main(argc - 1, argv + 1);
	}
	if(strcmp(argv[1], "target") == 0) {
		return target_main(argc - 1, argv + 1);
	}
	fprintf(stderr, "Unrecognized subcommand. See usage.\n");
	usage_exit(stderr, 1);
	return 0; // This never happens
}
