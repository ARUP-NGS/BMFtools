#include <unistd.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <limits.h>
#include <getopt.h>

#include "htslib/sam.h"
#include "samtools.h"
#include "khash.h"

HASH_MAP_INIT_INT64(fm, uint64_t)
HASH_MAP_INIT_INT64(rc, uint64_t)


typedef struct {
	uint64_t n_pass;
	uint64_t n_fail;
	uint64_t fm_sum;
	uint64_t rc_sum;
	khash_t(fm) *fm;
	khash_t(rc) *rc;
} famstats_t;

typedef struct famstat_settings {
	uint32_t minMQ;
	uint32_t minFM;
} famstat_settings_t;

static inline void famstat_loop(famstats_t *s, bam1_t *b, famstat_settings_t *settings)
{
	int FM = bam_aux2i(bam_aux_get(b, "FM"));
	if(b->core.qual < settings->minMQ || FM < settings->minFM) {
		++settings->n_fail;
		return;
	}
	++settings->n_pass;
	int RC = bam_aux2i(bam_aux_get(b, "RC"));
}


bamstats_t *famstat_core(samFile *fp, bam_hdr_t *h, famstat_settings_t *settings)
{
    famstat_t *s;
    bam1_t *b;
    bam1_core_t *c;
    int ret;
    s = (famstat_t*)calloc(1, sizeof(famstat_t));
    b = bam_init1();
    c = &b->core;
    while ((ret = sam_read1(fp, h, b)) >= 0)
        famstat_loop(s, c, settings);
    bam_destroy1(b);
    if (ret != -1)
        fprintf(stderr, "[famstat_core] Truncated file? Continue anyway.\n");
    return s;
}

static const char *percent(char *buffer, long long n, long long total)
{
    if (total != 0) sprintf(buffer, "%.2f%%", (float)n / total * 100.0);
    else strcpy(buffer, "N/A");
    return buffer;
}

static void usage_exit(FILE *fp, int exit_status)
{
    fprintf(fp, "Usage: samtools flagstat [--input-fmt-option OPT=VAL] <in.bam>\n");
    exit(exit_status);
}

int famstat(int argc, char *argv[])
{
    samFile *fp;
    bam_hdr_t *header;
    famstat_t *s;
    char b0[16], b1[16];
    hts_opt *in_opts = NULL;
    int c;

    enum {
        INPUT_FMT_OPTION = CHAR_MAX+1,
    };

    static const struct option lopts[] = {
        {"input-fmt-option",  required_argument, NULL, INPUT_FMT_OPTION},
        {NULL, 0, NULL, 0}
    };

    while ((c = getopt_long(argc, argv, "", lopts, NULL)) >= 0) {
        switch (c) {
        case INPUT_FMT_OPTION:
            if (hts_opt_add(&in_opts, optarg) < 0)
                usage_exit(stderr, EXIT_FAILURE);
            break;
        default:
            usage_exit(stderr, EXIT_FAILURE);
        }
    }

    if (argc != optind+1) {
        if (argc == optind) usage_exit(stdout, EXIT_SUCCESS);
        else usage_exit(stderr, EXIT_FAILURE);
    }
    fp = sam_open(argv[optind], "r");
    if (fp == NULL) {
        print_error_errno("flagstat", "Cannot open input file \"%s\"", argv[optind]);
        return 1;
    }
    if (hts_opt_apply(fp, in_opts)) {
        fprintf(stderr, "Failed to apply input-fmt-options\n");
        return 1;
    }

    if (hts_set_opt(fp, CRAM_OPT_REQUIRED_FIELDS,
                    SAM_FLAG | SAM_MAPQ | SAM_RNEXT)) {
        fprintf(stderr, "Failed to set CRAM_OPT_REQUIRED_FIELDS value\n");
        return 1;
    }

    if (hts_set_opt(fp, CRAM_OPT_DECODE_MD, 0)) {
        fprintf(stderr, "Failed to set CRAM_OPT_DECODE_MD value\n");
        return 1;
    }

    header = sam_hdr_read(fp);
    if (header == NULL) {
        fprintf(stderr, "Failed to read header for \"%s\"\n", argv[optind]);
        return 1;
    }
    s = famstat_core(fp, header);
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
    free(s);
    bam_hdr_destroy(header);
    sam_close(fp);
    hts_opt_free(in_opts);
    return 0;
}
