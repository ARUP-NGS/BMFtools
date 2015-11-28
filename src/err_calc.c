#include "err_calc.h"

void err_usage_exit(FILE *fp, int retcode)
{
	fprintf(fp, "Usage: Not written\n"
			"bmftools err -o <out.tsv> <reference.fasta> <input.csrt.bam>\n"
			"Opts:\n\t-h/-?\tThis helpful help menu!\n");
	exit(retcode);
}

// Add new dimension read 1/read 2 :'(
void make4d(FILE *fp, fullerr_t *e)
{
	for(uint32_t cycle = 0; cycle < e->l; ++cycle) {
		for(uint32_t qn = 0; qn < nqscores; ++qn) {
			fprintf(fp, "%i", e->r1->final[0][qn][cycle]);
			for(uint32_t bn = 1; bn < 4; ++bn)
				fprintf(fp, ":%i", e->r1->final[bn][qn][cycle]);
			if(qn != nqscores - 1) fprintf(fp, ",");
		}
		fprintf(fp, "|");
		for(uint32_t qn = 0; qn < nqscores; ++qn) {
			fprintf(fp, "%i", e->r2->final[0][qn][cycle]);
			for(uint32_t bn = 1; bn < 4; ++bn)
				fprintf(fp, ":%i", e->r2->final[bn][qn][cycle]);
			if(qn != nqscores - 1) fprintf(fp, ",");
		}
		fprintf(fp, "\n");
	}
}

void err_report(FILE *fp, fullerr_t *e)
{
	uint64_t min_obs = min_obs;
	fprintf(stderr, "Beginning error report.\n");
	fprintf(fp, "{\n{\"total_read\": %lu},\n{\"total_skipped\": %lu},\n", e->nread, e->nskipped);
	uint64_t n1_obs = 0, n1_err = 0, n1_ins = 0;
	uint64_t n2_obs = 0, n2_err = 0, n2_ins = 0;
	// n_ins is number with insufficient observations to report.
	for(int i = 0; i < 4; ++i) {
		for(int j = 0; j < nqscores; ++j) {
			for(int k = 0; k < e->l; ++k) {
				n1_obs += e->r1->obs[i][j][k];
				n2_obs += e->r2->obs[i][j][k];
				n1_err += e->r1->err[i][j][k];
				n2_err += e->r2->err[i][j][k];
				if(e->r1->rates[i][j][k] == DBL_MAX)
					++n1_ins;
				if(e->r2->rates[i][j][k] == DBL_MAX)
					++n2_ins;
			}
		}
	}
	uint64_t n_cases = nqscores * 4 * e->l;
	fprintf(stderr, "{\"read1\": {\"total_error\": %f},\n{\"total_obs\": %lu},\n{\"total_err\": %lu}"
			",\n{\"number_insufficient\": %lu},\n{\"n_cases\": %lu}},",
			(double)n1_err / n1_obs, n1_obs, n1_err, n1_ins, n_cases);
	fprintf(stderr, "{\"read2\": {\"total_error\": %f},\n{\"total_obs\": %lu},\n{\"total_err\": %lu}"
			",\n{\"number_insufficient\": %lu},\n{\"n_cases\": %lu}},",
			(double)n2_err / n2_obs, n2_obs, n2_err, n2_ins, n_cases);
	fprintf(fp, "}");
	return;
}

void readerr_destroy(readerr_t *e){
	for(int i = 0; i < 4; ++i) {
		for(int j = 0; j < nqscores; ++j) {
			free(e->obs[i][j]);
			free(e->err[i][j]);
			if(e->final[i][j]) free(e->final[i][j]);
			if(e->rates[i][j]) free(e->rates[i][j]);
		}
		free(e->obs[i]);
		free(e->err[i]);
		if(e->final[i]) free(e->final[i]);
		if(e->rates[i]) free(e->rates[i]);
		if(e->qrates[i]) free(e->qrates[i]);
		if(e->qcounts[i]) free(e->qcounts[i]);
	}
	if(e->qrates) free(e->qrates), e->qrates = NULL;
	if(e->qcounts) free(e->qcounts), e->qcounts = NULL;
	if(e->rates) free(e->rates), e->rates = NULL;
	free(e->obs);
	free(e->err);
	free(e);
}

void rate_calc(readerr_t *e)
{
	uint64_t min_obs = min_obs;
	e->rates = (double ***)malloc(4 * sizeof(double **));
	e->qrates = (double **)malloc(4 * sizeof(double *));
	e->qdiffs = (int **)malloc(4 * sizeof(int *));
	uint64_t **qobs = (uint64_t **)malloc(4 * sizeof(uint64_t *));
	uint64_t **qerr = (uint64_t **)malloc(4 * sizeof(uint64_t *));
	double **qsum = (double **)malloc(4 * sizeof(double *));
	uint64_t **qcounts = (uint64_t **)malloc(4 * sizeof(uint64_t *));
	for(int i = 0; i < 4; ++i) {
		qcounts[i] = (uint64_t *)calloc(e->l, sizeof(uint64_t));
		qsum[i] = (double *)calloc(e->l, sizeof(double));
		qobs[i] = (uint64_t *)calloc(e->l,  sizeof(uint64_t));
		qerr[i] = (uint64_t *)calloc(e->l,  sizeof(uint64_t));
		e->qrates[i] = (double *)malloc(e->l * sizeof(double));
		e->rates[i] = (double **)malloc(nqscores * sizeof(double *));
		for(int j = 0; j < nqscores; ++j) {
			e->rates[i][j] = (double *)malloc(e->l * sizeof(double));
			for(uint32_t l = 0; l < e->l; ++l) {
				e->rates[i][j][l] = (e->obs[i][j][l] >= min_obs) ? (double)e->err[i][j][l] / e->obs[i][j][l] : DBL_MAX;
				qsum[i][l] += (-10. * log10(j + 2)) * e->obs[i][j][l];
				qcounts += e->obs[i][j][l];
				qobs[i][l] += e->obs[i][j][l];
				qerr[i][l] += e->err[i][j][l];
			}
		}
	}
	int i;
	uint64_t l;
	for(i = 0; i < 4; ++i) {
		for(l = 0; l < e->l; ++l)
			qsum[i][l] /= qcounts[i][l];
		free(qcounts[i]);
	}
	free(qcounts);

	for(i = 0; i < 4; ++i) {
		for(l = 0; l < e->l; ++l) {
			e->qrates[i][l] = (qobs[i][l]) ? (double)qerr[i][l] / qobs[i][l]: DBL_MAX;
			e->qdiffs[i][l] = (e->qrates[i][l] != DBL_MAX) ?
					pv2ph(e->qrates[i][l]) - pv2ph(qsum[i][l]): 0.;
			// Difference between measured error rate and observed error rate
		}
		free(qobs[i]);
		free(qerr[i]);
	}
	e->final = (int ***)malloc(sizeof(int **) * 4);
	for(i = 0; i < 4; ++i) {
		e->final[i] = (int **)malloc(sizeof(int *) * nqscores);
		for(int j = 0; j < nqscores; ++j) {
			e->final[i][j] = (int *)malloc(sizeof(int) * e->l);
			for(l = 0; l < e->l; ++l)
				e->final[i][j][l] = (e->rates[i][j][l] != DBL_MAX) ?
						pv2ph(e->rates[i][j][l]): j + 2 + e->qdiffs[i][l];
		}
	}
	free(qobs);
	free(qerr);
	// Now impute with the diffs for those which had insufficient observations.
	return;
}
// TODO: finish cleanup. Delete the qrates and qdiff arrays



void err_core(char *fname, faidx_t *fai, fullerr_t *f, htsFormat *open_fmt)
{
	fprintf(stderr, "Opening bam file.\n");
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
		if(b->core.flag & 2816 || b->core.tid < 0) {// UNMAPPED, SECONDARY, SUPPLEMENTARY, QCFAIL
			++f->nskipped;
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

		if(++f->nread % 1000000 == 0)
			fprintf(stderr, "[%s] Records read: %lu.\n", __func__, f->nread);
#if !NDEBUG
		assert(b->core.tid >= 0);
#endif
		if(b->core.tid != last_tid) {
			if(ref) free(ref);
			fprintf(stderr, "Loading ref sequence for refernece contig with name %s.\n", hdr->target_name[b->core.tid]);
			ref = fai_fetch(fai, hdr->target_name[b->core.tid], &len);
			fprintf(stderr, "Finished loading ref sequence for contig '%s'.\n", hdr->target_name[b->core.tid]);
			last_tid = b->core.tid;
		}
		// rc -> read count
		// fc -> reference base count
		int i, ind, rc, fc;
		const readerr_t *r = (b->core.flag & BAM_FREAD1) ? f->r1: f->r2;
		fprintf(stderr, "Pointer to readerr_t r: %p.\n", r);
		const int32_t pos = b->core.pos;
		for(i = 0, rc = 0, fc = 0; i < b->core.n_cigar; ++i) {
			//fprintf(stderr, "Qual %p, seq %p, cigar %p.\n", seq, qual, cigar);
			uint8_t s;
			const uint32_t op = cigar[i];
			const uint32_t len = bam_cigar_oplen(op);
			switch(bam_cigar_op(op)) {
			case BAM_CMATCH:
				for(ind = 0; ind < len; ++ind) {
					s = bam_seqi(seq, ind + rc);
					if(s == HTS_N || ref[pos + fc + ind] == 'N') continue;
					if(qual[ind + rc] > nqscores + 1) { // nqscores + 2 - 1
						fprintf(stderr, "Quality score is too high. int: %i. char: %c. Max permitted: %lu.\n", (int)qual[ind + rc], qual[ind + rc], nqscores + 1);
						exit(EXIT_FAILURE);
					}
					if(bamseq2i[s] < 0) {
						fprintf(stderr, "HEY THIS SHOULDN'T BE NEGATIVE\n");
						exit(EXIT_FAILURE);
					}
					fprintf(stderr, "Incrementing obs %p, %p, %p, %lu\n", r->obs, r->obs[bamseq2i[s]], r->obs[bamseq2i[s]][qual[ind + rc] - 2],
							r->obs[bamseq2i[s]][qual[ind + rc] - 2][ind + rc]);
					fprintf(stderr, "Incrementing err %p, %p,\n", r->err, r->err[bamseq2i[s]]);
					++r->obs[bamseq2i[s]][qual[ind + rc] - 2][ind + rc];
					if(seq_nt16_table[(int)ref[pos + fc + ind]] != s) {
						//fprintf(stderr, "Found a mismatch before I died.\n");
						fprintf(stderr, "Incrementing err %p, %p, %p, %lu\n", r->err, r->err[bamseq2i[s]], r->err[bamseq2i[s]][qual[ind + rc] - 2],
								r->err[bamseq2i[s]][qual[ind + rc] - 2][ind + rc]);
							++r->err[bamseq2i[s]][qual[ind + rc] - 2][ind + rc];
					}
					//fprintf(stderr, "Finished incrementing.\n");
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
			// Default: break
			}
		}
	}
	fprintf(stderr, "Cleaning up after gathering my error data.\n");
	if(ref) free(ref);
	fprintf(stderr, "Deleting bam record.\n");
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

	FILE *ofp = NULL, *d2 = NULL, *d3 = NULL;
	while ((c = getopt(argc, argv, "m:o:h?")) >= 0) {
		switch (c) {
		case 'o': strcpy(outpath, optarg); break;
		case 'm': min_obs = atoi(optarg); break;
		case '3': d3 = fopen(optarg, "w"); break;
		case '2': d2 = fopen(optarg, "w"); break;
		case '?':
		case 'h':
			err_usage_exit(stderr, EXIT_SUCCESS);
		default:
			err_usage_exit(stderr, EXIT_FAILURE);
		}
	}


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
	fullerr_t *f = fullerr_init((size_t)b->core.l_qseq);
	sam_close(fp);
	fp = NULL;
	bam_destroy1(b);
	bam_hdr_destroy(header);
	header = NULL;
	err_core(argv[optind + 1], fai, f, &open_fmt);
	rate_calc(f->r1);
	rate_calc(f->r2);
	make4d(ofp, f);
	//err_report(ofp, e);
	fullerr_destroy(f);
	fclose(ofp);
	if(d3)
		fclose(d3);
	if(d2)
		fclose(d2);
	return 0;
}
