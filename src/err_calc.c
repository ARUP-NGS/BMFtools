#include "err_calc.h"

static uint64_t min_obs = 10;

void err_usage_exit(FILE *fp, int retcode)
{
	fprintf(fp, "Usage: Not written\n"
			"bmftools err -o <out.tsv> <reference.fasta> <input.csrt.bam>\n"
			"Opts:\n\t-h/-?\tThis helpful help menu!\n"
			"\t-3: Path to write the 3d offset array in tabular format.\n"
			"\t-4: Path to write the 4d measured array in tabular format.\n"
			);
	exit(retcode);
}

// Add new dimension read 1/read 2 :'(
void makefinal(FILE *fp, fullerr_t *e)
{
	fprintf(stderr, "Hey i'm making 4d!\n");
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
	fprintf(stderr, "Finishing makefinal.\n");
	return;
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
				if(e->r1->obs[i][j][k] < min_obs)
					++n1_ins;
				if(e->r2->obs[i][j][k] < min_obs)
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
			fprintf(stderr, "Destroying obs %i %i\n", i, j);
			if(e->obs[i][j]) free(e->obs[i][j]), e->obs[i][j] = NULL;
			fprintf(stderr, "Destroying err %i %i\n", i, j);
			if(e->err[i][j]) free(e->err[i][j]), e->err[i][j] = NULL;
			fprintf(stderr, "Destroying final %i, %i\n", i, j);
			if(e->final[i][j]) free(e->final[i][j]), e->final[i][j] = NULL;
			//if(e->rates[i][j]) free(e->rates[i][j]);
		}
		fprintf(stderr, "Destroying obs %i\n", i);
		if(e->obs[i]) free(e->obs[i]), e->obs[i] = NULL;
		fprintf(stderr, "Destroying err %i\n", i);
		if(e->err[i]) free(e->err[i]), e->err[i] = NULL;
		fprintf(stderr, "Destroying qobs %i\n", i);
		if(e->qobs[i]) free(e->qobs[i]), e->qobs[i] = NULL;
		fprintf(stderr, "Destroying qerr %i\n", i);
		if(e->qerr[i]) free(e->qerr[i]), e->qerr[i] = NULL;
		fprintf(stderr, "Destroying final %i\n", i);
		if(e->final[i]) free(e->final[i]), e->final[i] = NULL;
		fprintf(stderr, "Destroying qpvsum %i\n", i);
		if(e->qpvsum[i]) free(e->qpvsum[i]), e->qpvsum[i] = NULL;
		fprintf(stderr, "Destroying qdiffs %i\n", i);
		if(e->qdiffs[i]) free(e->qdiffs[i]), e->qdiffs[i] = NULL;
		//if(e->qrates && e->qrates[i]) free(e->qrates[i]);
		//if(e->qcounts && e->qcounts[i]) free(e->qcounts[i]);
	}
	//if(e->qrates) free(e->qrates), e->qrates = NULL;
	//if(e->qcounts) free(e->qcounts), e->qcounts = NULL;
	//if(e->rates) free(e->rates), e->rates = NULL;
	if(e->obs) free(e->obs), e->obs = NULL;
	if(e->err) free(e->err), e->err = NULL;
	if(e->qobs) free(e->qobs), e->qobs = NULL;
	if(e->qerr) free(e->qerr), e->qerr = NULL;
	if(e->final) free(e->final), e->final = NULL;
	if(e->qpvsum) free(e->qpvsum), e->qpvsum = NULL;
	if(e->qdiffs) free(e->qdiffs), e->qdiffs = NULL;
	free(e), e = NULL;
}


void err_core(char *fname, faidx_t *fai, fullerr_t *f, htsFormat *open_fmt)
{
	if(!f->r1)
		f->r1 = readerr_init(f->l);
	if(!f->r2)
		f->r2 = readerr_init(f->l);//
	samFile *fp = sam_open_format(fname, "r", open_fmt);
	bam_hdr_t *hdr = sam_hdr_read(fp);
	if (hdr == NULL) {
		exit(EXIT_FAILURE);
	}
	int r, len;
	int32_t last_tid = -1;
	bam1_t *b = bam_init1();
	char *ref = NULL; // Will hold the sequence for a chromosome
	while((r = sam_read1(fp, hdr, b)) != -1) {
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
		static const int bamseq2i[] = {-1, 0, 1, -1, 2, -1, -1, -1, 3};
		if(b->core.tid != last_tid) {
			if(ref) free(ref);
			fprintf(stderr, "Loading ref sequence for contig with name %s.\n", hdr->target_name[b->core.tid]);
			ref = fai_fetch(fai, hdr->target_name[b->core.tid], &len);
			fprintf(stderr, "Finished loading ref sequence for contig '%s'.\n", hdr->target_name[b->core.tid]);
			last_tid = b->core.tid;
		}
		// rc -> read count
		// fc -> reference base count
		int i, ind, rc, fc;
		const readerr_t *r = (b->core.flag & BAM_FREAD1) ? f->r1: f->r2;
		//fprintf(stderr, "Pointer to readerr_t r: %p.\n", r);
		const int32_t pos = b->core.pos;
		for(i = 0, rc = 0, fc = 0; i < b->core.n_cigar; ++i) {
			//fprintf(stderr, "Qual %p, seq %p, cigar %p.\n", seq, qual, cigar);
			int s; // seq value, base index
			const uint32_t op = cigar[i];
			const uint32_t len = bam_cigar_oplen(op);
			switch(bam_cigar_op(op)) {
			case BAM_CMATCH:
				for(ind = 0; ind < len; ++ind) {
					s = bam_seqi(seq, ind + rc);
					//fprintf(stderr, "Bi value: %i. s: %i.\n", bi, s);
					if(s == HTS_N || ref[pos + fc + ind] == 'N') continue;
					if(qual[ind + rc] > nqscores + 1) { // nqscores + 2 - 1
						fprintf(stderr, "Quality score is too high. int: %i. char: %c. Max permitted: %lu.\n", (int)qual[ind + rc], qual[ind + rc], nqscores + 1);
						exit(EXIT_FAILURE);
					}
					/*
					fprintf(stderr, "indices: %i, ind + rec %i, qual val %i\n", bi, ind + rc, qual[ind + rc] - 2);
					fprintf(stderr, "Incrementing obs %p, %p, %i\n", r->obs, r->obs[bi], qual[ind + rc] - 2);
					fprintf(stderr, "Incrementing err %p, %p, %i\n", r->err, r->err[bi], qual[ind + rc] - 2);
					fprintf(stderr, "Incrementing obs %p, %p, %p\n", r->obs, r->obs[bi], r->obs[bi][qual[ind + rc] - 2]);
					fprintf(stderr, "Incrementing err %p, %p, %p\n", r->err, r->err[bi], r->err[bi][qual[ind + rc] - 2]);
					fprintf(stderr, "Incrementing obs %p, %p, %p, %lu\n", r->obs, r->obs[bi], r->obs[bi][qual[ind + rc] - 2], r->obs[bi][qual[ind + rc] - 2][ind + rc]);
					fprintf(stderr, "Incrementing err %p, %p, %p, %lu\n", r->err, r->err[bi], r->err[bi][qual[ind + rc] - 2], r->err[bi][qual[ind + rc] - 2][ind + rc]);
					*/
					++r->obs[bamseq2i[s]][qual[ind + rc] - 2][ind + rc];
					if(seq_nt16_table[(int)ref[pos + fc + ind]] != s)
						//fprintf(stderr, "Found a mismatch before I died.\n");
						//fprintf(stderr, "Incrementing err %p, %p, %p, %lu\n", r->err, r->err[bamseq2i[s]], r->err[bamseq2i[s]][qual[ind + rc] - 2],
						//		r->err[bamseq2i[s]][qual[ind + rc] - 2][ind + rc]);
							++r->err[bamseq2i[s]][qual[ind + rc] - 2][ind + rc];
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
#if !NDEBUG
	fprintf(stderr, "Deleting bam record at pointer %p.\n", b);
#else
	fprintf(stderr, "Deleting bam record.\n");
#endif
	bam_destroy1(b);
	return;
}

void write_4d_errs(FILE *fp, fullerr_t *f)
{
	/*
	fprintf(fp, "##Each line is a cycle on the sequencer"
			"#R1/A/2\tR1/C/2\tR1/G/2\t&c.\t\n");
	for(uint64_t k = 0; k < f->l; ++k) {
		for(int j = 0; j < nqscores; ++j) {
			for(int i = 0; i < 4; ++i) {
				fprintf(fp, (i) ? ":%0.8f": "%0.8f", f->r1->rates[i][j][k]);
			}
			if(j != nqscores - 1) fprintf(fp, ",");
		}
		fprintf(fp, "|");
		for(int j = 0; j < nqscores; ++j) {
			for(int i = 0; i < 4; ++i) {
				fprintf(fp, (i) ? ":%0.8f": "%0.8f", f->r2->rates[i][j][k]);
			}
			if(j != nqscores - 1) fprintf(fp, ",");
		}
		fprintf(fp, "\n");
	}
	*/
}

void fill_qvals(fullerr_t *f)
{
	fprintf(stderr, "Beginning fill_qvals\n");
	for(int i = 0; i < 4; ++i) {
		for(uint64_t l = 0; l < f->l; ++l) {
			for(int j = 1; j < nqscores; ++j) { // Skip qualities of 2
				f->r1->qpvsum[i][l] +=  pow(10., (double)(-0.1 * (j + 2))) * f->r1->obs[i][j][l];
				f->r2->qpvsum[i][l] +=  pow(10., (double)(-0.1 * (j + 2))) * f->r2->obs[i][j][l];
				f->r1->qobs[i][l] += f->r1->obs[i][j][l];
				f->r2->qobs[i][l] += f->r2->obs[i][j][l];
				f->r1->qerr[i][l] += f->r1->err[i][j][l];
				f->r2->qerr[i][l] += f->r2->err[i][j][l];
			}
		}
	}
	for(int i = 0; i < 4; ++i) {
		for(uint64_t l = 0; l < f->l; ++l) {
			f->r1->qpvsum[i][l] /= f->r1->qobs[i][l]; // Get average ILMN-reported quality
			f->r2->qpvsum[i][l] /= f->r2->qobs[i][l]; // Divide by observations of cycle/base call
			f->r1->qdiffs[i][l] = pv2ph((double)f->r1->qerr[i][l] / f->r1->qobs[i][l]) - pv2ph(f->r1->qpvsum[i][l]);
			f->r2->qdiffs[i][l] = pv2ph((double)f->r2->qerr[i][l] / f->r2->qobs[i][l]) - pv2ph(f->r2->qpvsum[i][l]);
			if(f->r1->qobs[i][l] < min_obs) f->r1->qdiffs[i][l] = 0;
			if(f->r2->qobs[i][l] < min_obs) f->r2->qdiffs[i][l] = 0;
		}
	}
	FILE *qvh = fopen("qvalues.txt", "w");
	for(uint64_t l = 0; l < f->l; ++l) {
		for(int i = 0; i < 4; ++i){
			fprintf(qvh, i ? ":%i": "%i", f->r1->qdiffs[i][l]);
			fputc('|', qvh);
			fprintf(qvh, i ? ":%i": "%i", f->r2->qdiffs[i][l]);
		}
		fputc('\n', qvh);
	}
	fclose(qvh);
}

void write_counts(fullerr_t *f, FILE *cp, FILE *ep)
{
	const char *bstr = "ACGT";
	FILE *dictwrite = fopen("dict.txt", "w");
	fprintf(dictwrite, "{\n\t");
	for(uint32_t l = 0; l < f->l; ++l) {
		for(int j = 0; j < nqscores; ++j) {
			for(int i = 0; i < 4; ++i) {
				fprintf(dictwrite, "'r1,%c,%i,%u,obs': %lu,\n\t", bstr[i], j + 2, l + 1, f->r1->obs[i][j][l]);
				fprintf(dictwrite, "'r2,%c,%i,%u,obs': %lu,\n\t", bstr[i], j + 2, l + 1, f->r2->obs[i][j][l]);
				fprintf(dictwrite, "'r1,%c,%i,%u,err': %lu,\n\t", bstr[i], j + 2, l + 1, f->r1->err[i][j][l]);
				if(i == 3 && j == nqscores - 1 && l == f->l - 1)
					fprintf(dictwrite, "'r2,%c,%i,%u,err': %lu\n}", bstr[i], j + 2, l + 1, f->r2->err[i][j][l]);
				else
					fprintf(dictwrite, "'r2,%c,%i,%u,err': %lu,\n\t", bstr[i], j + 2, l + 1, f->r2->err[i][j][l]);
				fprintf(cp, i ? ":%lu": "%lu", f->r1->obs[i][j][l]);
				fprintf(ep, i ? ":%lu": "%lu", f->r1->err[i][j][l]);
			}
			if(j != nqscores - 1)
				fprintf(ep, ","), fprintf(cp, ",");
		}
		fprintf(ep, "|"), fprintf(cp, "|");
		for(int j = 0; j < nqscores; ++j) {
			for(int i = 0; i < 4; ++i) {
				fprintf(cp, i ? ":%lu": "%lu", f->r2->obs[i][j][l]);
				fprintf(ep, i ? ":%lu": "%lu", f->r2->err[i][j][l]);
			}
			if(j != nqscores - 1)
				fprintf(ep, ","), fprintf(cp, ",");
		}
		fprintf(ep, "\n"), fprintf(cp, "\n");
	}
	fclose(dictwrite);
	return;
}

void write_3d_offsets(FILE *fp, fullerr_t *f)
{
	for(uint64_t k = 0; k < f->l; ++k) {
		for(int i = 0; i < 4; ++i) {
			fprintf(fp, i ? ":%i": "%i", f->r1->qdiffs[i][k]);
		}
		fprintf(fp, "|");
		for(int i = 0; i < 4; ++i) {
			fprintf(fp, i ? ":%i": "%i", f->r2->qdiffs[i][k]);
		}
		fprintf(fp, "\n");
	}
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

	FILE *ofp = NULL, *d4 = NULL, *d3 = NULL;
	while ((c = getopt(argc, argv, "4:3:m:o:h?")) >= 0) {
		switch (c) {
		case 'o': strcpy(outpath, optarg); break;
		case 'm': min_obs = atoi(optarg); break;
		case '4': d4 = fopen(optarg, "w"); break;
		case '3': d3 = fopen(optarg, "w"); break;
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
	fprintf(stderr, "Core finished.\n");
	FILE *ch = fopen("counts.txt", "w"),*eh = fopen("errs.txt", "w");
	write_counts(f, ch, eh);
	fill_qvals(f);
	cfclose(ch); cfclose(eh);
	/*
	rate_calc(f->r1);
	rate_calc(f->r2);
	makefinal(ofp, f);
	*/
	if(d4)
		fprintf(stderr, "Writin' 4d errs\n"),
		write_4d_errs(d4, f), fclose(d4);
	if(d3)
		fprintf(stderr, "Writin' 3d errs\n"),
		write_3d_offsets(d3, f), fclose(d3);
	//err_report(ofp, e);
	fullerr_destroy(f);
	fclose(ofp);
	return 0;
}
