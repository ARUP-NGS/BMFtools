#include <unistd.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <limits.h>
#include <getopt.h>
#include <inttypes.h>

#include "htslib/sam.h"
//#include "samtools.h"
#include "khash.h"

KHASH_MAP_INIT_INT64(fm, uint64_t)
KHASH_MAP_INIT_INT64(rc, uint64_t)


typedef struct famstats {
	uint64_t n_pass;
	uint64_t n_fail;
	uint64_t allfm_sum;
	uint64_t allfm_counts;
	uint64_t allrc_sum;
	uint64_t allrc_counts;
	uint64_t realfm_sum;
	uint64_t realfm_counts;
	uint64_t realrc_sum;
	uint64_t realrc_counts;
	khash_t(fm) *fm;
	khash_t(rc) *rc;
	khiter_t ki;
	int khr;
	uint8_t *data;
} famstats_t;

typedef struct famstat_settings {
	uint32_t minMQ;
	uint32_t minFM;
} famstat_settings_t;

static inline void print_stats(famstats_t *stats)
{
	fprintf(stderr, "Number passing filters: %"PRIu64".\n", stats->n_pass);
	fprintf(stderr, "Number failing filters: %"PRIu64".\n", stats->n_fail);
	fprintf(stderr, "Summed FM (total founding reads): %"PRIu64".\n", stats->allfm_sum);
	fprintf(stderr, "Summed RC(total reverse-complemented reads): %"PRIu64".\n", stats->allrc_sum);
}

static inline void tag_test(uint8_t *data, const char *tag)
{
	if(data)
		return;
	fprintf(stderr, "Required bam tag '%s' not found. Abort mission!\n", tag);
	exit(EXIT_FAILURE);
}


static inline void famstat_loop(famstats_t *s, bam1_t *b, famstat_settings_t *settings)
{
	++s->n_pass;
	uint8_t *data;
	data = bam_aux_get(b, "FM");
	tag_test(data, "FM");
	int FM = bam_aux2i(data);
#if !NDEBUG
	fprintf(stderr, "FM tag: %i.\n", FM);
#endif
	if(b->core.qual < settings->minMQ || FM < settings->minFM) {
		++s->n_fail;
		return;
	}
	data = bam_aux_get(b, "RC");
	tag_test(data, "RC");
	int RC = bam_aux2i(data);
#if !NDEBUG
	fprintf(stderr, "RC tag: %i.\n", RC);
#endif
	if(FM > 1) {
		++s->realfm_counts; s->realfm_sum += FM; ++s->realrc_counts; s->realrc_sum += RC;
	}
	++s->allfm_counts; s->allfm_sum += FM; ++s->allrc_counts; s->allrc_sum += RC;
	s->ki = kh_get(fm, s->fm, FM);
	if(s->ki == kh_end(s->fm)) {
		s->ki = kh_put(fm, s->fm, FM, &s->khr);
		kh_val(s->fm, s->ki) = 1;
	}
	else {
		++kh_val(s->fm, s->ki);
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
	while ((ret = sam_read1(fp, h, b)) >= 0) {
		famstat_loop(s, b, settings);
		if(!(++count % 250000))
			fprintf(stderr, "Number of records processed: %"PRIu64".\n", count);
	}
	bam_destroy1(b);
	if (ret != -1)
		fprintf(stderr, "[famstat_core] Truncated file? Continue anyway.\n");
	return s;
}

static void usage_exit(FILE *fp, int exit_status)
{
	fprintf(fp, "Usage: famstat <in.bam>\n");
	fprintf(fp, "Opts:\n-m Set minimum mapping quality. Default: 0.\n");
	fprintf(fp, "-f Set minimum family size. Default: 0.\n");
	exit(exit_status);
}

int main(int argc, char *argv[])
{
	samFile *fp;
	bam_hdr_t *header;
	famstats_t *s;
	int c;

	famstat_settings_t *settings = (famstat_settings_t *)calloc(1, sizeof(famstat_settings_t));
	settings->minMQ = 0;
	settings->minFM = 0;

	while ((c = getopt(argc, argv, "m:f:")) >= 0) {
		switch (c) {
		case 'm':
			settings->minMQ = atoi(optarg); break;
			break;
		case 'f':
			settings->minFM = atoi(optarg); break;
			break;
		case 'h':
			usage_exit(stderr, EXIT_SUCCESS);
		default:
			usage_exit(stderr, EXIT_FAILURE);
		}
	}
	fprintf(stderr, "[fam_stat]: Running main with minMQ %i and minFM %i.\n", settings->minMQ, settings->minFM);

	if (argc != optind+1) {
		if (argc == optind) usage_exit(stdout, EXIT_SUCCESS);
		else usage_exit(stderr, EXIT_FAILURE);
	}
	fp = sam_open(argv[optind], "r");
	if (fp == NULL) {
		fprintf(stderr, "Cannot open input file \"%s\"", argv[optind]);
		exit(EXIT_FAILURE);
	}

	header = sam_hdr_read(fp);
	if (header == NULL) {
		fprintf(stderr, "Failed to read header for \"%s\"\n", argv[optind]);
		exit(EXIT_FAILURE);
	}
	s = famstat_core(fp, header, settings);
	/*
	printf("%lld + %lld in total (QC-passed reads + QC-failed reads)\n", s->n_reads[0], s->n_reads[1]);
	printf("%lld + %lld secondary\n", s->n_secondary[0], s->n_secondary[1]);
	printf("%lld + %lld supplementary\n", s->n_supp[0], s->n_supp[1]);
	printf("%lld + %lld duplicates\n", s->n_dup[0], s->n_dup[1]);
	printf("%lld + %lld mapped (%s : %s)\n", s->n_mapped[0], s->n_mapped[1], percent(b0, s->n_mapped[0], s->n_reads[0]), percent(b1, s->n_mapped[1], s->n_reads[1]));
	printf("%lld + %lld paired in sequencing\n", s->n_pair_all[0], s->n_pair_all[1]);
	printf("%lld + %lld read1\n", s->n_read1[0], s->n_read1[1]);
	printf("%lld + %lld read2\n", s->n_read2[0], s->n_read2[1]);
	printf("%lld + %lld properly paired (%s : %s)\n", s->n_pair_good[0], s->n_pair_good[1], percent(b0, s->n_pair_good[0], s->n_pair_all[0]), percent(b1, s->n_pair_good[1], s->n_pair_all[1]));
	printf("%lld + %lld with itself and mate mapped\n", s->n_pair_map[0], s->n_pair_map[1]);
	printf("%lld + %lld singletons (%s : %s)\n", s->n_sgltn[0], s->n_sgltn[1], percent(b0, s->n_sgltn[0], s->n_pair_all[0]), percent(b1, s->n_sgltn[1], s->n_pair_all[1]));
	printf("%lld + %lld with mate mapped to a different chr\n", s->n_diffchr[0], s->n_diffchr[1]);
	printf("%lld + %lld with mate mapped to a different chr (mapQ>=5)\n", s->n_diffhigh[0], s->n_diffhigh[1]);
	*/
	double mean_fm_all = (double)(s->allfm_sum) / s->allfm_counts;
	double mean_fm_real = (double)(s->realfm_sum) / s->realfm_counts;
	double mean_rcfrac_all = (double)(s->allfm_sum) / s->allrc_sum;
	double mean_rcfrac_real = (double)(s->realfm_sum) / s->realrc_sum;
	fprintf(stdout, "Mean Family Size (all)\t%f\n", mean_fm_all);
	fprintf(stdout, "Mean Family Size (FM > 1)\t%f\n", mean_fm_real);
	fprintf(stdout, "Mean Reverse Complement Fraction (all)\t%f.\n", mean_rcfrac_all);
	fprintf(stdout, "Mean Reverse Complement Fraction (FM > 1)\t%f\n", mean_rcfrac_real);
	fprintf(stdout, "Number of families total\t%"PRIu64"\n", s->allfm_counts);
	fprintf(stdout, "Number of families (FM > 1)\t%"PRIu64"\n", s->realfm_counts);
	fprintf(stdout, "Number of raw reads reverse-complemented\t%"PRIu64"\n", s->allrc_counts);
	fprintf(stdout, "Number of raw reads reverse-complemented in families (FM > 1)\t%"PRIu64"\n", s->realrc_counts);
	free(s);
	free(settings);
	bam_hdr_destroy(header);
	sam_close(fp);
	return 0;
}
