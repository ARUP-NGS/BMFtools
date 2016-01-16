#include "bmf_famstats.h"

int RVWarn = 1;

int frac_usage_exit(FILE *fp, int code) {
	fprintf(fp, "bmftools famstats frac <opts> <in.bam>\n"
			"Opts:\n-m minFM to accept. REQUIRED.\n"
			"-h, -?: Return usage.\n");
	exit(code);
	return code; // This never happens
}

typedef struct {
	uint64_t n; // Number of times observed
	uint64_t fm; // Number of times observed
} fm_t;

int fm_cmp(const void *fm1, const void *fm2)
{
	return ((fm_t *)fm1)->fm - ((fm_t *)fm2)->fm;
}


static void print_hashstats(famstats_t *stats, FILE *fp)
{
	fm_t *fms = (fm_t *)malloc(sizeof(fm_t) * MAX2(stats->fm->n_occupied, stats->rc->n_occupied));
	int i = 0;
	fprintf(fp, "#Family size\tNumber of families\n");
	for(stats->ki = kh_begin(stats->fm); stats->ki != kh_end(stats->fm); ++stats->ki) {
		if(!kh_exist(stats->fm, stats->ki)) continue;
		fms[i++] = (fm_t){.fm = kh_key(stats->fm, stats->ki), .n = kh_val(stats->fm, stats->ki)};
	}
	qsort(fms, stats->fm->n_occupied, sizeof(fm_t), fm_cmp);
	for(i = 0; i < stats->fm->n_occupied; ++i) fprintf(fp, "%lu\t%lu\n", fms[i].fm, fms[i].n);
	fprintf(fp, "#RV'd in family\tNumber of families\n");
	for(stats->ki = kh_begin(stats->rc); stats->ki != kh_end(stats->rc); ++stats->ki) {
		if(!kh_exist(stats->rc, stats->ki)) continue;
		fms[i++] = (fm_t){.fm = kh_key(stats->rc, stats->ki), .n = kh_val(stats->rc, stats->ki)};
	}
	qsort(fms, stats->rc->n_occupied, sizeof(fm_t), fm_cmp);
	for(i = 0; i < stats->rc->n_occupied; ++i) fprintf(fp, "%lu\t%lu\n", fms[i].fm, fms[i].n);
	free(fms);
}


int target_usage_exit(FILE *fp, int success)
{
	fprintf(fp, "Usage: bmftools famstats target <opts> <in.bam>\nOpts:\n-b Path to bed file.\n"
			"-p padding. Number of bases around bed regions to pad. Default: 25.\n"
			"-h, -?: Return usage.\n");
	exit(success);
	return success;
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

	if(argc < 4) return target_usage_exit(stderr, EXIT_SUCCESS);

	if(strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0) return target_usage_exit(stderr, EXIT_SUCCESS);

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
			return target_usage_exit(stderr, EXIT_SUCCESS);
		}
	}


	if(padding == (uint32_t)-1) {
		padding = 25;
		LOG_INFO("Padding not set. Set to 25.\n");
	}

	if (argc != optind+1)
		return (argc == optind) ? target_usage_exit(stdout, EXIT_SUCCESS): target_usage_exit(stderr, EXIT_FAILURE);

	if(!bedpath) {
		fprintf(stderr, "[E:%s] Bed path required for famstats target. See usage.\n", __func__);
		return target_usage_exit(stderr, EXIT_FAILURE);
	}

	fp = sam_open(argv[optind], "r");
	if (fp == NULL) {
		LOG_ERROR("Cannot open input file \"%s\"", argv[optind]);
	}

	header = sam_hdr_read(fp);
	if (!header) {
		LOG_ERROR("Failed to read header for \"%s\"\n", argv[optind]);
	}
	bed = parse_bed_hash(bedpath, header, padding);
	uint64_t fm_target = 0, total_fm = 0, count = 0;
	bam1_t *b = bam_init1();
	while (LIKELY((c = sam_read1(fp, header, b)) >= 0)) {
		target_loop(b, bed, &fm_target, &total_fm);
		if(UNLIKELY(++count % notification_interval == 0)) fprintf(stderr, "[%s] Number of records processed: %lu.\n", __func__, count);
	}
	bam_destroy1(b);
	bam_hdr_destroy(header);
	sam_close(fp);
	bed_destroy_hash(bed);
	free(bedpath);
	fprintf(stdout, "Fraction of raw reads on target: %f. \nTotal raw reads: %lu. Raw reads on target: %lu.\n",
			(double)fm_target / total_fm, total_fm, fm_target);
	return EXIT_SUCCESS;
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
	if(stats->dr_counts) {
		fprintf(fp, "#Duplex fraction of unique observations\t%0.12lf\n", (double)stats->dr_counts / stats->n_pass);
		fprintf(fp, "#Fraction of raw reads in duplex families\t%0.12lf\n", (double)stats->dr_sum / stats->allfm_sum);
		fprintf(fp, "#Mean fraction of reverse reads within each duplex family\t%0.12lf\n", stats->dr_rc_frac_sum / stats->dr_rc_sum);
		fprintf(fp, "#Mean fraction of reverse reads within all duplex families\t%0.12lf\n", (double)stats->dr_rc_sum / stats->dr_sum);
	}
	print_hashstats(stats, fp);
}

static inline void tag_test(const uint8_t *data, const char *tag)
{
	if(UNLIKELY(!data))
		fprintf(stderr, "[E:%s] Required bam tag '%s' not found. Abort mission!\n", __func__, tag),
		exit(EXIT_FAILURE);
}


static inline void famstat_loop(famstats_t *s, bam1_t *b, famstat_settings_t *settings)
{
	//tag_test(data, "FM");
	if((b->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY | BAM_FQCFAIL | BAM_FREAD2)) ||
			b->core.qual < settings->minMQ) return;
	const int FM = bam_aux2i(bam_aux_get(b, "FM"));
	const int RV = bam_aux2i(bam_aux_get(b, "RV"));
	if(FM < settings->minFM || !bam_aux2i(bam_aux_get(b, "FP"))) {
		++s->n_fail;
		return;
	}
	++s->n_pass;

	if(FM > 1) ++s->realfm_counts, s->realfm_sum += FM, s->realrc_sum += RV;
	++s->allfm_counts, s->allfm_sum += FM, s->allrc_sum += RV;

	LOG_DEBUG("RV value: %i. allrc_sum: %lu. realrc_sum: %lu.\n", RV, s->allrc_sum, s->realrc_sum);

	if((s->ki = kh_get(fm, s->fm, FM)) == kh_end(s->fm))
		s->ki = kh_put(fm, s->fm, FM, &s->khr), kh_val(s->fm, s->ki) = 1;
	else ++kh_val(s->fm, s->ki);
	if((s->ki = kh_get(rc, s->rc, RV)) == kh_end(s->rc))
		s->ki = kh_put(rc, s->rc, RV, &s->khr), kh_val(s->rc, s->ki) = 1;
	else ++kh_val(s->rc, s->ki);

	uint8_t *dr_data = bam_aux_get(b, "DR");
	if(dr_data && bam_aux2i(dr_data)) {
		s->dr_sum += FM;
		++s->dr_counts;
		s->dr_rc_sum += RV;
		s->dr_rc_frac_sum += (double)RV / FM;
	}
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
	ret = sam_read1(fp, h, b);
	if(ret < 0) {
		LOG_ERROR("Could not read from input bam %s. Abort!\n", fp->fn);
	}
	++count;
	famstat_loop(s, b, settings);
	check_bam_tag(b, "FP");
	check_bam_tag(b, "RV");
	check_bam_tag(b, "FM");
	while (LIKELY((ret = sam_read1(fp, h, b)) >= 0)) {
		famstat_loop(s, b, settings);
		if(UNLIKELY(++count % settings->notification_interval == 0))
			LOG_INFO("Number of records processed: %lu.\n", count);
	}
	bam_destroy1(b);
	if (ret != -1) LOG_WARNING("Truncated file? Continue anyway.\n");
	return s;
}

static int usage_exit(FILE *fp, int exit_status)
{
	fprintf(fp, "Usage: bmftools famstats\n");
	fprintf(fp, "Subcommands: \nfm\tFamily Size stats\n"
			"frac\tFraction of raw reads in family sizes >= minFM parameter.\n"
			"target\tFraction of raw reads on target.\n");
	exit(exit_status);
	return exit_status;
}

static int fm_usage_exit(FILE *fp, int exit_status)
{
	fprintf(fp, "Usage: bmftools famstats fm <opts> <in.bam>\n");
	fprintf(fp, "-m Set minimum mapping quality. Default: 0.\n");
	fprintf(fp, "-f Set minimum family size. Default: 0.\n");
	exit(exit_status);
	return exit_status;
}

int fm_main(int argc, char *argv[])
{
	samFile *fp;
	bam_hdr_t *header;
	famstats_t *s;
	int c;


	famstat_settings_t *settings = (famstat_settings_t *)calloc(1, sizeof(famstat_settings_t));
	settings->notification_interval = 1000000uL;

	while ((c = getopt(argc, argv, "m:f:n:h?")) >= 0) {
		switch (c) {
		case 'm':
			settings->minMQ = atoi(optarg); break;
			break;
		case 'f':
			settings->minFM = atoi(optarg); break;
			break;
		case 'n': settings->notification_interval = strtoull(optarg, NULL, 0); break;
		case '?': // Fall-through!
		case 'h':
			return fm_usage_exit(stderr, EXIT_SUCCESS);
		}
	}

	if (argc != optind+1) {
		if (argc == optind) fm_usage_exit(stdout, EXIT_SUCCESS);
		else fm_usage_exit(stderr, EXIT_FAILURE);
	}

	LOG_INFO("Running main with minMQ %i and minFM %i.\n", settings->minMQ, settings->minFM);

	fp = sam_open(argv[optind], "r");
	if (!fp) {
		LOG_ERROR("Cannot open input file \"%s\"", argv[optind]);
	}

	header = sam_hdr_read(fp);
	if (!header) {
		LOG_ERROR("Failed to read header for \"%s\"\n", argv[optind]);
	}
	s = famstat_core(fp, header, settings);
	print_stats(s, stdout);
	kh_destroy(fm, s->fm);
	kh_destroy(rc, s->rc);
	free(s);
	free(settings);
	bam_hdr_destroy(header);
	sam_close(fp);
	return EXIT_SUCCESS;
}

static inline void frac_loop(bam1_t *b, int minFM, uint64_t *fm_above, uint64_t *fm_total)
{
	if((b->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY | BAM_FQCFAIL | BAM_FREAD2)) ||
			bam_aux2i(bam_aux_get(b, "FP")) == 0)
		return;
	const int FM = bam_aux2i(bam_aux_get(b, "FM"));
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
		LOG_ERROR("minFM not set. frac_main meaningless without it. Result: 1.0.\n");
	}
	LOG_INFO("Running frac main minFM %i.\n", minFM);

	if (argc != optind+1) {
		if (argc == optind) frac_usage_exit(stdout, EXIT_SUCCESS);
		else frac_usage_exit(stderr, EXIT_FAILURE);
	}
	fp = sam_open(argv[optind], "r");
	if (!fp) {
		LOG_ERROR("Cannot open input file \"%s\".\n", argv[optind]);
	}

	header = sam_hdr_read(fp);
	if (!header) {
		LOG_ERROR("Failed to read header for \"%s\".\n", argv[optind]);
	}
	uint64_t fm_above = 0, total_fm = 0, count = 0;
	bam1_t *b = bam_init1();
	// Check to see if the required tags are present before starting.
	if((c = sam_read1(fp, header, b)) < 0) {
		LOG_ERROR("Could not read initial record from input file '%s'. Abort!\n", argv[optind]);
	}
	check_bam_tag(b, "FP");
	check_bam_tag(b, "RV");
	check_bam_tag(b, "FM");
	check_bam_tag(b, "FA");
	frac_loop(b, minFM, &fm_above, &total_fm), ++count;
	while (LIKELY(c = sam_read1(fp, header, b)) >= 0) {
		frac_loop(b, minFM, &fm_above, &total_fm);
		if(UNLIKELY(!(++count % notification_interval)))
			fprintf(stderr, "[%s] Number of records processed: %lu.\n", __func__, count);
	}
	fprintf(stdout, "#Fraction of raw reads with >= minFM %i: %f.\n", minFM, (double)fm_above / total_fm);
	bam_destroy1(b);
	bam_hdr_destroy(header);
	sam_close(fp);
	return EXIT_SUCCESS;
}

int famstats_main(int argc, char *argv[])
{
	if(argc < 2)
		return usage_exit(stderr, EXIT_FAILURE);
	if(strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0)
		return usage_exit(stderr, EXIT_SUCCESS);
	if(strcmp(argv[1], "fm") == 0)
		return fm_main(argc - 1, argv + 1);
	if(strcmp(argv[1], "frac") == 0)
		return frac_main(argc - 1, argv + 1);
	if(strcmp(argv[1], "target") == 0)
		return target_main(argc - 1, argv + 1);
	fprintf(stderr, "[E:%s] Unrecognized subcommand '%s'. See usage.\n", __func__, argv[1]);
	return usage_exit(stderr, EXIT_FAILURE);
	return EXIT_FAILURE; // This never happens
}
