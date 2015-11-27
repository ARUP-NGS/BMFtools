#include "err_calc.h"

int min_obs = 200;

void err_usage_exit(FILE *fp, int retcode)
{
	fprintf(fp, "Usage: Not written\n"
			"bmftools err -o <out.tsv> <reference.fasta> <input.csrt.bam>\n"
			"Opts:\n\t-h/-?\tThis helpful help menu!\n");
	exit(retcode);
}

void err_report(FILE *fp, errcnt_t *e)
{
	fprintf(fp, "{\n{\"total_read\": %lu},\n{\"total_skipped\": %lu},\n", e->nread, e->nskipped);
	uint64_t n_obs = 0, n_err = 0;
	for(int i = 0; i < 4; ++i) {
		for(int j = 0; j < nqscores; ++j) {
			for(int k = 0; k < e->l; ++k) {
				n_obs += e->obs[i][j][k];
				n_err += e->err[i][j][k];
				if(e->obs[i][j][k] >= min_obs) {
					e->rates[i][j][k] = (double)e->err[i][j][k] / e->obs[i][j][k];
				}
			}
		}
	}
	fprintf(fp, "{\"total_error\": %lf},\n{\"total_obs\": %lu},\n{\"total_err\": %lu}",
			(double)n_err / n_obs, n_obs, n_err);
	fprintf(fp, "}");
	return;
}

static int bamseq2i(uint8_t seqi) {
	fprintf(stderr, "Now calling bamseq2i.\n");
	switch(seqi) {
	case HTS_A:
		return 0;
	case HTS_C:
		return 1;
	case HTS_G:
		return 2;
	default: // HTS_T, since HTS_N was already checked for
		return 3;
	}
}

void err_core(char *fname, faidx_t *fai, errcnt_t *e, htsFormat *open_fmt)
{
	samFile *fp = sam_open_format(fname, "r", open_fmt);
	bam_hdr_t *hdr = sam_hdr_read(fp);
	if (hdr == NULL) {
		fprintf(stderr, "[famstat_err_main]: Failed to read header for \"%s\"\n", fname);
		exit(EXIT_FAILURE);
	}
	int r, len;
	int32_t last_tid = -1;
	bam1_t *b = bam_init1();
	char *ref = NULL; // Will hold the sequence for a chromosome
	while((r = sam_read1(fp, hdr, b)) != -1) {
		fprintf(stderr, "Read new record with name %s.\n", bam_get_qname(b));
		if(b->core.flag & 2816 || b->core.tid < 0) { // UNMAPPED, SECONDARY, SUPPLEMENTARY, QCFAIL
			++e->nskipped;
			continue;
		}
		const uint8_t *seq = (uint8_t *)bam_get_seq(b);
		const uint8_t *qual = (uint8_t *)bam_get_qual(b);
		const uint32_t *cigar = bam_get_cigar(b);
		if(!cigar) {
			fprintf(stderr, "Could not get bam cigar. Abort!\n");
			exit(EXIT_FAILURE);
		}
		if(!seq) {
			fprintf(stderr, "Could not get bam seq. Abort!\n");
			exit(EXIT_FAILURE);
		}
		if(!qual) {
			fprintf(stderr, "Could not get bam qual. Abort!\n");
			exit(EXIT_FAILURE);
		}

		if(++e->nread % 10 == 0)
			fprintf(stderr, "[%s] Records read: %lu.\n", __func__, e->nread);
#if !NDEBUG
		assert(b->core.tid >= 0);
#endif
		if(b->core.tid != last_tid) {
			if(ref) free(ref);
			fprintf(stderr, "Loading ref sequence for tid %i, name %s.\n", b->core.tid, hdr->target_name[b->core.tid]);
			ref = fai_fetch(fai, hdr->target_name[b->core.tid], &len);
			fprintf(stderr, "Finished loading ref sequence for tid %i, name %s.\n", b->core.tid, hdr->target_name[b->core.tid]);
			last_tid = b->core.tid;
		}
		// rc -> read count
		// fc -> reference base count
		int i, ind, rc, fc;
		const int32_t pos = b->core.pos;
		fprintf(stderr, "Calculation time.\n");
		for(i = 0, rc = 0, fc = 0; i < b->core.n_cigar; ++i) {
			fprintf(stderr, "Qual %p, seq %p, cigar %p.\n", seq, qual, cigar);
			uint8_t s;
			const uint32_t op = cigar[i];
			const uint32_t len = bam_cigar_oplen(op);
			fprintf(stderr, "Got to the switch!\n");
			switch(bam_cigar_op(op)) {
			case BAM_CMATCH:
				for(ind = 0; ind < len; ++ind) {
					s = bam_seqi(seq, ind + rc);
					if(s == HTS_N) continue;
					if(qual[ind + rc] > nqscores - 1) {
						fprintf(stderr, "Quality score is too high. int: %i. char: %c. Max permitted: %i.\n", (int)qual[ind + rc], qual[ind + rc], nqscores - 1);
						exit(EXIT_FAILURE);
					}
#if !NDEBUG
					if(ind + rc > )
					fprintf(stderr, "Pointer: qual ind %i, qual val %i, %p.\n", ind + rc, qual[ind + rc] - 2, e->obs[bamseq2i(s)]);
					fprintf(stderr, "Incrementing observations.\n");
#endif
					++e->obs[bamseq2i(s)][qual[ind + rc] - 2][ind + rc];
					if(seq_nt16_table[(int)ref[pos + fc + ind]] != s) {
						fprintf(stderr, "Found a mismatch before I died.\n");
						++e->err[bamseq2i(s)][qual[ind + rc] - 2][ind + rc];
					}
					fprintf(stderr, "Finished incrementing.\n");
				}
			case BAM_CEQUAL:
			case BAM_CDIFF:
				rc += len;
				fc += len;
				break;
			case BAM_CSOFT_CLIP:
			case BAM_CHARD_CLIP:
			case BAM_CINS:
				rc += len;
				break;
			case BAM_CREF_SKIP:
			case BAM_CDEL:
				fc += len;
				break;
#if !NDEBUG
			default:
				fprintf(stderr, "Default case thrown for op %i", bam_cigar_op(op)); break;
#endif
			// Default: break
			}
#if !NDEBUG
			fprintf(stderr, "Keep going! Finished cigar operation #%i (1-based).\n", i + 1); break;
#endif
		}
		fprintf(stderr, "Added things to the thingamajigger.\n");
	}
	if(ref) free(ref);
	bam_destroy1(b);
	return;
}

int err_main(int argc, char *argv[])
{
	htsFormat open_fmt;
	memset(&open_fmt, 0, sizeof(htsFormat));
	open_fmt.category = sequence_data;
	open_fmt.format = bam;
	open_fmt.version.major = 1;
	open_fmt.version.minor = 3;
	samFile *fp = NULL;
	bam_hdr_t *header = NULL;
	int c;
	char outpath[500] = "";

	if(argc < 2) {
		err_usage_exit(stderr, EXIT_FAILURE);
	}

	if(strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0) err_usage_exit(stderr, EXIT_SUCCESS);

	while ((c = getopt(argc, argv, "o:h?")) >= 0) {
		switch (c) {
		case 'o': strcpy(outpath, optarg); break;
		case '?':
		case 'h':
			err_usage_exit(stderr, EXIT_SUCCESS);
		default:
			err_usage_exit(stderr, EXIT_FAILURE);
		}
	}

	FILE *ofp = NULL;
	ofp = (outpath[0]) ? fopen(outpath, "w"): stdout;

	if (argc != optind+2)
		err_usage_exit(stderr, EXIT_FAILURE);

	faidx_t *fai = fai_load(argv[optind]);

	fp = sam_open_format(argv[optind + 1], "r", &open_fmt);
	if (fp == NULL) {
		fprintf(stderr, "[famstat_err_main]: Cannot open input file \"%s\"", argv[optind]);
		exit(EXIT_FAILURE);
	}

	header = sam_hdr_read(fp);
	if (header == NULL) {
		fprintf(stderr, "[famstat_err_main]: Failed to read header for \"%s\"\n", argv[optind]);
		exit(EXIT_FAILURE);
	}
#if !NDEBUG
	//for(int i = 0; i < header->n_targets; ++i)
		//fprintf(stderr, "Target name %i: %s\n", i, header->target_name[i]);
#endif
	// Get read length from the first read
	bam1_t *b = bam_init1();
	c = sam_read1(fp, header, b);
	errcnt_t *e = errcnt_init((size_t)b->core.l_qseq);
	sam_close(fp);
	fp = NULL;
	bam_destroy1(b);
	err_core(argv[optind + 1], fai, e, &open_fmt);
	err_report(ofp, e);
	bam_destroy1(b);
	bam_hdr_destroy(header);
	sam_close(fp);
	errcnt_destroy(e);
	fclose(ofp);
	return 0;
}
