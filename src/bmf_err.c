#include "bmf_err.h"

#define min_obs 1000uL

void err_usage_exit(FILE *fp, int retcode)
{
	fprintf(fp,
			"Usage: bmftools err -o <out.tsv> <reference.fasta> <input.csrt.bam>\n"
			"Flags:\n"
			"-h/-?\t\tThis helpful help menu!\n"
			"-r:\t\tName of contig. If set, only reads aligned to this contig are considered\n"
			"-3:\t\tPath to write the 3d offset array in tabular format.\n"
			"-f:\t\tPath to write the full measured error rates in tabular format.\n"
			"-b:\t\tPath to write the cycle/base call error rates in tabular format.\n"
			"-c:\t\tPath to write the cycle error rates in tabular format.\n"
			);
	exit(retcode);
}

void write_final(FILE *fp, fullerr_t *e)
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
}

void readerr_destroy(readerr_t *e){
	for(int i = 0; i < 4; ++i) {
		for(int j = 0; j < nqscores; ++j) {
			cond_free(e->obs[i][j]);
			cond_free(e->err[i][j]);
			cond_free(e->final[i][j]);
		}
		cond_free(e->obs[i]);
		cond_free(e->err[i]);
		cond_free(e->qobs[i]);
		cond_free(e->qerr[i]);
		cond_free(e->final[i]);
		cond_free(e->qpvsum[i]);
		cond_free(e->qdiffs[i]);
	}
	cond_free(e->obs);
	cond_free(e->err);
	cond_free(e->qerr);
	cond_free(e->qobs);
	cond_free(e->final);
	cond_free(e->qpvsum);
	cond_free(e->qdiffs);
	cond_free(e);
}


void err_core(char *fname, faidx_t *fai, fullerr_t *f, htsFormat *open_fmt)
{
	if(!f->r1) f->r1 = readerr_init(f->l);
	if(!f->r2) f->r2 = readerr_init(f->l);
	samFile *fp = sam_open_format(fname, "r", open_fmt);
	bam_hdr_t *hdr = sam_hdr_read(fp);
	if (!hdr) {
		fprintf(stderr, "[E:%s] Failed to read input header from bam %s. Abort!\n", __func__, fname);
		exit(EXIT_FAILURE);
	}
	int r, len;
	int32_t last_tid = -1;
	bam1_t *b = bam_init1();
	char *ref = NULL; // Will hold the sequence for a chromosome
	int tid_to_study = -1;
	if(f->refcontig) {
		for(int i = 0; i < hdr->n_targets; ++i) {
			if(!strcmp(hdr->target_name[i], f->refcontig)) {
				tid_to_study = i; break;
			}
		}
		if(tid_to_study < 0) {
			fprintf(stderr, "Contig %s not found in bam header. Abort mission!\n", f->refcontig);
			exit(EXIT_FAILURE);
		}
	}
	while(LIKELY((r = sam_read1(fp, hdr, b)) != -1)) {
		if((b->core.flag & 2820) || (f->refcontig && tid_to_study != b->core.tid)) {++f->nskipped; continue;} // UNMAPPED, SECONDARY, SUPPLEMENTARY, QCFAIL
		const uint8_t *seq = (uint8_t *)bam_get_seq(b);
		const uint8_t *qual = (uint8_t *)bam_get_qual(b);
		const uint32_t *cigar = bam_get_cigar(b);
#if !NDEBUG
		ifn_abort(cigar);
		ifn_abort(seq);
		ifn_abort(qual);
#endif

		if(++f->nread % 1000000 == 0) fprintf(stderr, "[%s] Records read: %lu.\n", __func__, f->nread);
#if !NDEBUG
		assert(b->core.tid >= 0);
#endif
		if(b->core.tid != last_tid) {
			cond_free(ref);
			fprintf(stderr, "[%s] Loading ref sequence for contig with name %s.\n", __func__, hdr->target_name[b->core.tid]);
			ref = fai_fetch(fai, hdr->target_name[b->core.tid], &len);
			last_tid = b->core.tid;
		}
		const readerr_t *r = (b->core.flag & BAM_FREAD1) ? f->r1: f->r2;
		const int32_t pos = b->core.pos;
		for(int i = 0, rc = 0, fc = 0; i < b->core.n_cigar; ++i) {
			int s; // seq value, base index
			const uint32_t len = bam_cigar_oplen(*cigar);
			switch(bam_cigar_op(*cigar++)) {
			case BAM_CMATCH:
			case BAM_CEQUAL:
			case BAM_CDIFF:
				for(int ind = 0; ind < len; ++ind) {
					s = bam_seqi(seq, ind + rc);
					//fprintf(stderr, "Bi value: %i. s: %i.\n", bi, s);
					if(s == HTS_N || ref[pos + fc + ind] == 'N') continue;
#if !NDEBUG
					if(UNLIKELY(qual[ind + rc] > nqscores + 1)) { // nqscores + 2 - 1
						fprintf(stderr, "[E:%s] Quality score is too high. int: %i. char: %c. Max permitted: %lu.\n",
								__func__, (int)qual[ind + rc], qual[ind + rc], nqscores + 1);
						exit(EXIT_FAILURE);
					}
#endif
					++r->obs[bamseq2i[s]][qual[ind + rc] - 2][ind + rc];
					if(seq_nt16_table[(int8_t)ref[pos + fc + ind]] != s) ++r->err[bamseq2i[s]][qual[ind + rc] - 2][ind + rc];
				}
				rc += len; fc += len;
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
			}
		}
	}
	fprintf(stderr, "[D:%s] Cleaning up after gathering my error data.\n", __func__);
	cond_free(ref);
	bam_destroy1(b);
	bam_hdr_destroy(hdr), sam_close(fp);
}

void err_core_se(char *fname, faidx_t *fai, fullerr_t *f, htsFormat *open_fmt)
{
	if(!f->r1) f->r1 = readerr_init(f->l);
	samFile *fp = sam_open_format(fname, "r", open_fmt);
	bam_hdr_t *hdr = sam_hdr_read(fp);
	if (!hdr) {
		fprintf(stderr, "[E:%s] Failed to read input header from bam %s. Abort!\n", __func__, fname);
		exit(EXIT_FAILURE);
	}
	int len;
	int32_t last_tid = -1;
	bam1_t *b = bam_init1();
	char *ref = NULL; // Will hold the sequence for a chromosome
	int tid_to_study = -1;
	const readerr_t *rerr = f->r1;
	if(f->refcontig) {
		for(int i = 0; i < hdr->n_targets; ++i) {
			if(!strcmp(hdr->target_name[i], f->refcontig)) {
				tid_to_study = i; break;
			}
		}
		if(tid_to_study < 0) {
			fprintf(stderr, "Contig %s not found in bam header. Abort mission!\n", f->refcontig);
			exit(EXIT_FAILURE);
		}
	}
	int c;
	while(LIKELY((c = sam_read1(fp, hdr, b)) != -1)) {
		if((b->core.flag & 2820) || (f->refcontig && tid_to_study != b->core.tid)) {++f->nskipped; continue;} // UNMAPPED, SECONDARY, SUPPLEMENTARY, QCFAIL
		const uint8_t *seq = (uint8_t *)bam_get_seq(b);
		const uint8_t *qual = (uint8_t *)bam_get_qual(b);
		const uint32_t *cigar = bam_get_cigar(b);
#if !NDEBUG
		ifn_abort(cigar);
		ifn_abort(seq);
		ifn_abort(qual);
#endif

		if(++f->nread % 1000000 == 0) fprintf(stderr, "[%s] Records read: %lu.\n", __func__, f->nread);
#if !NDEBUG
		assert(b->core.tid >= 0);
#endif
		if(b->core.tid != last_tid) {
			cond_free(ref);
			fprintf(stderr, "[%s] Loading ref sequence for contig with name %s.\n", __func__, hdr->target_name[b->core.tid]);
			ref = fai_fetch(fai, hdr->target_name[b->core.tid], &len);
			last_tid = b->core.tid;
		}
		// rc -> read count
		// fc -> reference base count
		//fprintf(stderr, "Pointer to readerr_t r: %p.\n", r);
		const int32_t pos = b->core.pos;
		for(int i = 0, rc = 0, fc = 0; i < b->core.n_cigar; ++i) {
			//fprintf(stderr, "Qual %p, seq %p, cigar %p.\n", seq, qual, cigar);
			int s; // seq value, base index
			const uint32_t len = bam_cigar_oplen(*cigar);
			switch(bam_cigar_op(*cigar++)) {
			case BAM_CMATCH:
			case BAM_CEQUAL:
			case BAM_CDIFF:
				for(int ind = 0; ind < len; ++ind) {
					s = bam_seqi(seq, ind + rc);
					//fprintf(stderr, "Bi value: %i. s: %i.\n", bi, s);
					if(s == HTS_N || ref[pos + fc + ind] == 'N') continue;
#if !NDEBUG
					if(UNLIKELY(qual[ind + rc] > nqscores + 1)) { // nqscores + 2 - 1
						fprintf(stderr, "[E:%s] Quality score is too high. int: %i. char: %c. Max permitted: %lu.\n",
								__func__, (int)qual[ind + rc], qual[ind + rc], nqscores + 1);
						exit(EXIT_FAILURE);
					}
#endif
					++rerr->obs[bamseq2i[s]][qual[ind + rc] - 2][ind + rc];
					if(seq_nt16_table[(int8_t)ref[pos + fc + ind]] != s) ++rerr->err[bamseq2i[s]][qual[ind + rc] - 2][ind + rc];
				}
				rc += len; fc += len;
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
			}
		}
	}
	fprintf(stderr, "[D:%s] Cleaning up after gathering my error data.\n", __func__);
	cond_free(ref);
	bam_destroy1(b);
	bam_hdr_destroy(hdr), sam_close(fp);
}


void write_full_rates(FILE *fp, fullerr_t *f)
{
	uint64_t l;
	int i, j;
	for(l = 0; l < f->l; ++l) {
		for(j = 0; j < nqscores; ++j) {
			for(i = 0; i < 4; ++i) {
				if(f->r1->obs[i][j][l])
					fprintf(fp, i ? ":%0.12f": "%0.12f", (double)f->r1->err[i][j][l] / f->r1->obs[i][j][l]);
				else fputs(i ? ":-1337": "-1337", fp);
			}
			if(j != nqscores - 1) fputc(',', fp);
		}
		fputc('|', fp);
		for(j = 0; j < nqscores; ++j) {
			for(i = 0; i < 4; ++i) {
				if(f->r2->obs[i][j][l])
					fprintf(fp, i ? ":%0.12f": "%0.12f", (double)f->r2->err[i][j][l] / f->r2->obs[i][j][l]);
				else
					fputs(i ? ":-1337": "-1337", fp);
			}
			if(j != nqscores - 1) fputc(',', fp);
		}
		fputc('\n', fp);
	}
}



void write_base_rates(FILE *fp, fullerr_t *f)
{
	fputs("#Cycle\tR1A\tR1C\tR1G\tR1T\tR2A\tR2C\tR2G\tR2T\n", fp);
	for(uint64_t l = 0; l < f->l; ++l) {
		int i;
		fprintf(fp, "%lu\t", l + 1);
		for(i = 0; i < 4; ++i)
			fprintf(fp, i ? "\t%0.12f": "%0.12f", (double)f->r1->qobs[i][l] / f->r1->qerr[i][l]);
		fputc('|', fp);
		for(i = 0; i < 4; ++i)
			fprintf(fp, i ? "\t%0.12f": "%0.12f", (double)f->r2->qobs[i][l] / f->r2->qerr[i][l]);
		fputc('\n', fp);
	}
}

void write_cycle_rates(FILE *fp, fullerr_t *f)
{
	fputs("#Cycle\tRead 1 Error Rate\tRead 2 Error Rate\n", fp);
	for(uint64_t l = 0; l < f->l; ++l) {
		fprintf(fp, "%lu\t", l + 1);
		int sum1 = 0, sum2 = 0, counts1 = 0, counts2 = 0;
		for(int i = 0; i < 4; ++i) {
			sum1 += f->r1->qerr[i][l];
			counts1 += f->r1->qerr[i][l];
			sum2 += f->r2->qerr[i][l];
			counts2 += f->r2->qerr[i][l];
		}
		fprintf(fp, "%0.12f\t", (double)sum1 / counts1);
		fprintf(fp, "%0.12f\n", (double)sum2 / counts2);
	}
}

void impute_scores(fullerr_t *f)
{
	int j, i;
	uint64_t l;
	for(i = 0; i < 4; ++i)
		for(l = 0; l < f->l; ++l)
			for(j = 0; j < nqscores; ++j)
				f->r1->final[i][j][l] = f->r1->qdiffs[i][l] + j + 2 > 0 ? f->r1->qdiffs[i][l] + j + 2: 0,
				f->r2->final[i][j][l] = f->r2->qdiffs[i][l] + j + 2 > 0 ? f->r2->qdiffs[i][l] + j + 2: 0;
}

void fill_qvals(fullerr_t *f)
{
	int i;
	uint64_t l;
	for(i = 0; i < 4; ++i) {
		for(l = 0; l < f->l; ++l) {
			for(int j = 1; j < nqscores; ++j) { // Skip qualities of 2
				f->r1->qpvsum[i][l] +=  pow(10., (double)(-0.1 * (j + 2))) * f->r1->obs[i][j][l];
				f->r2->qpvsum[i][l] +=  pow(10., (double)(-0.1 * (j + 2))) * f->r2->obs[i][j][l];
				f->r1->qobs[i][l] += f->r1->obs[i][j][l]; f->r2->qobs[i][l] += f->r2->obs[i][j][l];
				f->r1->qerr[i][l] += f->r1->err[i][j][l]; f->r2->qerr[i][l] += f->r2->err[i][j][l];
			}
		}
	}
	for(i = 0; i < 4; ++i) {
		for(l = 0; l < f->l; ++l) {
			f->r1->qpvsum[i][l] /= f->r1->qobs[i][l]; // Get average ILMN-reported quality
			f->r2->qpvsum[i][l] /= f->r2->qobs[i][l]; // Divide by observations of cycle/base call
			f->r1->qdiffs[i][l] = pv2ph((double)f->r1->qerr[i][l] / f->r1->qobs[i][l]) - pv2ph(f->r1->qpvsum[i][l]);
			f->r2->qdiffs[i][l] = pv2ph((double)f->r2->qerr[i][l] / f->r2->qobs[i][l]) - pv2ph(f->r2->qpvsum[i][l]);
			//fprintf(stderr, "qdiffs i, l (measured) is R1:%i R2:%i.\n", f->r1->qdiffs[i][l], f->r2->qdiffs[i][l]);
			if(f->r1->qobs[i][l] < min_obs) f->r1->qdiffs[i][l] = 0;
			if(f->r2->qobs[i][l] < min_obs) f->r2->qdiffs[i][l] = 0;
			//fprintf(stderr, "qdiffs %i, %lu after checking for %lu %lu > %lu min_obs is R1:%i R2:%i.\n", i, l, f->r1->qobs[i][l], f->r2->qobs[i][l], min_obs, f->r1->qdiffs[i][l], f->r2->qdiffs[i][l]);
		}
	}
}

void fill_sufficient_obs(fullerr_t *f)
{
#if DELETE_ME
	FILE *before_fs = fopen("before_fill_sufficient.txt", "w");
	write_3d_offsets(before_fs, f);
	fclose(before_fs);
#endif
	for(int i = 0; i < 4; ++i) {
		for(int j = 0; j < nqscores; ++j) {
			for(uint64_t l = 0; l < f->l; ++l) {
				if(f->r1->obs[i][j][l] > min_obs)
					f->r1->final[i][j][l] = pv2ph((double)f->r1->err[i][j][l] / f->r1->obs[i][j][l]);
				if(f->r2->obs[i][j][l] > min_obs)
					f->r2->final[i][j][l] = pv2ph((double)f->r2->err[i][j][l] / f->r2->obs[i][j][l]);
			}
		}
	}
#if DELETE_ME
	FILE *after_fs = fopen("after_fill_sufficient.txt", "w");
	write_3d_offsets(after_fs, f);
	fclose(after_fs);
#endif
}

void write_counts(fullerr_t *f, FILE *cp, FILE *ep)
{
	const char *const bstr = "ACGT";
	FILE *dictwrite = fopen("dict.txt", "w");
	fprintf(dictwrite, "{\n\t");
	int i, j;
	uint32_t l;
	for(l = 0; l < f->l; ++l) {
		for(j = 0; j < nqscores; ++j) {
			for(i = 0; i < 4; ++i) {
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
		for(j = 0; j < nqscores; ++j) {
			for(i = 0; i < 4; ++i) {
				fprintf(cp, i ? ":%lu": "%lu", f->r2->obs[i][j][l]);
				fprintf(ep, i ? ":%lu": "%lu", f->r2->err[i][j][l]);
			}
			if(j != nqscores - 1)
				fprintf(ep, ","), fprintf(cp, ",");
		}
		fprintf(ep, "\n"), fprintf(cp, "\n");
	}
	fclose(dictwrite);
}

void write_3d_offsets(FILE *fp, fullerr_t *f)
{
	fprintf(fp, "#Cycle\tR1A\tR1C\tR1G\tR1T\tR2A\tR2C\tR2G\tR2T\n");
	for(uint64_t l = 0; l < f->l; ++l) {
		fprintf(fp, "%lu\t", l + 1);
		int i;
		for(i = 0; i < 4; ++i) fprintf(fp, i ? "\t%i": "%i", f->r1->qdiffs[i][l]);
		fputc('|', fp);
		for(i = 0; i < 4; ++i) fprintf(fp, i ? "\t%i": "%i", f->r2->qdiffs[i][l]);
		fputc('\n', fp);
	}
	return;
}

readerr_t *readerr_init(size_t l) {
	readerr_t *ret = (readerr_t *)calloc(1, sizeof(readerr_t));
	arr3d_init(ret->obs, l, uint64_t);
	arr3d_init(ret->err, l, uint64_t);
	arr3d_init(ret->final, l, int);
	arr2d_init(ret->qdiffs, l, int);
	arr2d_init(ret->qpvsum, l, double);
	arr2d_init(ret->qobs, l, uint64_t);
	arr2d_init(ret->qerr, l, uint64_t);
	ret->l = l;
	return ret;
}

fullerr_t *fullerr_init(size_t l) {
	fullerr_t *ret = (fullerr_t *)calloc(1, sizeof(fullerr_t));
	ret->l = l;
	ret->r1 = readerr_init(l);
	ret->r2 = readerr_init(l);
	return ret;
}

void fullerr_destroy(fullerr_t *e) {
	if(e->r1) readerr_destroy(e->r1), e->r1 = NULL;
	if(e->r2) readerr_destroy(e->r2), e->r2 = NULL;
	if(e->refcontig) free(e->refcontig), e->refcontig = NULL;
	free(e);
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

	FILE *ofp = NULL, *d3 = NULL, *df = NULL, *db = NULL, *dc = NULL;
	char refcontig[200] = "";
	while ((c = getopt(argc, argv, "r:c:b:f:3:o:h?")) >= 0) {
		switch (c) {
		case 'f': df = fopen(optarg, "w"); break;
		case 'o': strcpy(outpath, optarg); break;
		case '3': d3 = fopen(optarg, "w"); break;
		case 'c': dc = fopen(optarg, "w"); break;
		case 'b': db = fopen(optarg, "w"); break;
		case 'r': strcpy(refcontig, optarg); break;
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
	if(*refcontig) f->refcontig = strdup(refcontig);
	bam_hdr_destroy(header);
	header = NULL;
	err_core(argv[optind + 1], fai, f, &open_fmt);
	fprintf(stderr, "Core finished.\n");
	fai_destroy(fai);
	fill_qvals(f);
	impute_scores(f);
	fill_sufficient_obs(f);
	write_final(ofp, f);
	if(d3)
		fprintf(stderr, "Writin' 3d offsets.\n"),
		write_3d_offsets(d3, f), fclose(d3), d3 = NULL;
	if(df)
		fprintf(stderr, "Writin' read/base call/qscore/cycle error rates.\n"),
		write_full_rates(df, f), fclose(df), df = NULL;
	if(db)
		fprintf(stderr, "Writin' read/base call/cycle error rates.\n"),
		write_base_rates(db, f), fclose(db), db = NULL;
	if(dc)
		fprintf(stderr, "Writin' cycle error rates.\n"),
		write_cycle_rates(dc, f), fclose(dc), dc = NULL;
	fullerr_destroy(f);
	fclose(ofp);
	return 0;
}
