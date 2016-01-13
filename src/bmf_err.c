#include "bmf_err.h"

#define min_obs 1000uL

int err_usage(FILE *fp, int retcode)
{
	fprintf(fp,
			"Usage: bmftools err -o <out.tsv> <reference.fasta> <input.csrt.bam>\n"
			"Flags:\n"
			"-h/-?\t\tThis helpful help menu!\n"
			"-r:\t\tName of contig. If set, only reads aligned to this contig are considered\n"
			"-3:\t\tPath to write the 3d offset array in tabular format.\n"
			"-f:\t\tPath to write the full measured error rates in tabular format.\n"
			"-n:\t\tPath to write the cycle/nucleotide call error rates in tabular format.\n"
			"-c:\t\tPath to write the cycle error rates in tabular format.\n"
			"-b:\t\tPath to bed file for restricting analysis.\n"
			"-m:\t\tMinimum family size for inclusion. Default: 0.\n"
			"-M:\t\tMaximum family size for inclusion. Default: %i.\n",
			INT_MAX);
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

void fm_err_core(char *fname, faidx_t *fai, fmerr_t *f, htsFormat *open_fmt)
{
	samFile *fp = sam_open_format(fname, "r", open_fmt);
	bam_hdr_t *hdr = sam_hdr_read(fp);
	if (!hdr) {
		LOG_ERROR("Failed to read input header from bam %s. Abort!\n", fname);
	}
	int r, len, khr;
	int32_t last_tid = -1;
	bam1_t *b = bam_init1();
	char *ref = NULL; // Will hold the sequence for a chromosome
	int tid_to_study = -1;
	khiter_t k;
	if(f->refcontig) {
		for(int i = 0; i < hdr->n_targets; ++i) {
			if(!strcmp(hdr->target_name[i], f->refcontig)) {
				tid_to_study = i; break;
			}
		}
		if(tid_to_study < 0) {
			LOG_ERROR("Contig %s not found in bam header. Abort mission!\n", f->refcontig);
		}
	}
	while(LIKELY((r = sam_read1(fp, hdr, b)) != -1)) {
		const int FM = bam_aux2i(bam_aux_get(b, "FM"));
		const int RV = bam_aux2i(bam_aux_get(b, "RV"));
		if((b->core.flag & 2820) || (f->refcontig && tid_to_study != b->core.tid) ||
			(f->bed && bed_test(b, f->bed) == 0) || // Outside of region
			((f->flag & REQUIRE_DUPLEX) && (RV == FM || RV == 0)) || // Requires duplex
			(bam_aux2i(bam_aux_get(b, "FP")) == 0) // Fails barcode QC
			) {++f->nskipped; continue;} // UNMAPPED, SECONDARY, SUPPLEMENTARY, QCFAIL
		const uint8_t *seq = (uint8_t *)bam_get_seq(b);
		const uint8_t *qual = (uint8_t *)bam_get_qual(b);
		const uint32_t *cigar = bam_get_cigar(b);
#if !NDEBUG
		ifn_abort(cigar);
		ifn_abort(seq);
		ifn_abort(qual);
#endif

		if(++f->nread % 1000000 == 0) fprintf(stderr, "[%s] Records read: %lu.\n", __func__, f->nread);
		if(b->core.tid != last_tid) {
			cond_free(ref);
			LOG_DEBUG("Loading ref sequence for contig with name %s.\n", hdr->target_name[b->core.tid]);
			ref = fai_fetch(fai, hdr->target_name[b->core.tid], &len);
			if(!ref) {
				LOG_ERROR("Failed to load ref sequence for contig '%s'. Abort!\n", hdr->target_name[b->core.tid]);
			}
			last_tid = b->core.tid;
		}
		const int32_t pos = b->core.pos;
		k = kh_get(obs, f->hash, FM);
		if(k == kh_end(f->hash)) {
			k = kh_put(obs, f->hash, FM, &khr);
		}
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
						LOG_ERROR("Quality score is too high. int: %i. char: %c. Max permitted: %lu.\n",
								(int)qual[ind + rc], qual[ind + rc], nqscores + 1);
					}
#endif
					++kh_val(f->hash, k).obs;
					if(seq_nt16_table[(int8_t)ref[pos + fc + ind]] != s) ++kh_val(f->hash, k).err;
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
	cond_free(ref);
	bam_destroy1(b);
	bam_hdr_destroy(hdr), sam_close(fp);
}


void err_core(char *fname, faidx_t *fai, fullerr_t *f, htsFormat *open_fmt)
{
	if(!f->r1) f->r1 = readerr_init(f->l);
	if(!f->r2) f->r2 = readerr_init(f->l);
	samFile *fp = sam_open_format(fname, "r", open_fmt);
	bam_hdr_t *hdr = sam_hdr_read(fp);
	if (!hdr) {
		LOG_ERROR("Failed to read input header from bam %s. Abort!\n", fname);
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
			LOG_ERROR("Contig %s not found in bam header. Abort mission!\n", f->refcontig);
		}
	}
	while(LIKELY((r = sam_read1(fp, hdr, b)) != -1)) {
		const int FM = bam_aux2i(bam_aux_get(b, "FM"));
		const int RV = bam_aux2i(bam_aux_get(b, "RV"));
		if((b->core.flag & 2820) || (f->refcontig && tid_to_study != b->core.tid) ||
			(f->bed && bed_test(b, f->bed) == 0) || // Outside of region
			(FM < f->minFM) || (FM > f->maxFM) || // minFM outside of range
			((f->flag & REQUIRE_DUPLEX) && (RV == FM || RV == 0)) || // Requires duplex
			(bam_aux2i(bam_aux_get(b, "FP")) == 0) // Fails barcode QC
			) {++f->nskipped; continue;} // UNMAPPED, SECONDARY, SUPPLEMENTARY, QCFAIL
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
			LOG_DEBUG("Loading ref sequence for contig with name %s.\n", hdr->target_name[b->core.tid]);
			ref = fai_fetch(fai, hdr->target_name[b->core.tid], &len);
			if(!ref) {
				LOG_ERROR("[Failed to load ref sequence for contig '%s'. Abort!\n", hdr->target_name[b->core.tid]);
			}
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
						LOG_ERROR("Quality score is too high. int: %i. char: %c. Max permitted: %lu.\n",
								(int)qual[ind + rc], qual[ind + rc], nqscores + 1);
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
			LOG_ERROR("Contig %s not found in bam header. Abort mission!\n", f->refcontig);
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

		if(UNLIKELY(++f->nread % 1000000 == 0)) fprintf(stderr, "[%s] Records read: %lu.\n", __func__, f->nread);
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
		for(i = 0; i < 4; ++i) {
#if !NDEBUG
			LOG_DEBUG("obs: %lu. err: %lu.\n", f->r1->qerr[i][l], f->r1->qobs[i][l]);
#endif
			fprintf(fp, i ? "\t%0.12f": "%0.12f", (double)f->r1->qerr[i][l] / f->r1->qobs[i][l]);
		}
		fputc('|', fp);
		for(i = 0; i < 4; ++i)
			fprintf(fp, i ? "\t%0.12f": "%0.12f", (double)f->r2->qerr[i][l] / f->r2->qobs[i][l]);
		fputc('\n', fp);
	}
}


void write_global_rates(FILE *fp, fullerr_t *f)
{
	fprintf(fp, "##Parameters: minFM %i. maxFM %i.", f->minFM, f->maxFM);
	fputs("Duplex required: ", fp);
	fputs((f->flag & REQUIRE_DUPLEX) ? "True": "False", fp);
	fputc('\n', fp);
	int sum1 = 0, sum2 = 0, counts1 = 0, counts2 = 0;
	for(uint64_t l = 0; l < f->l; ++l) {
		for(int i = 0; i < 4; ++i) {
			sum1 += f->r1->qerr[i][l];
			counts1 += f->r1->qobs[i][l];
			sum2 += f->r2->qerr[i][l];
			counts2 += f->r2->qobs[i][l];
		}
	}
	fprintf(fp, "#Global Error Rates\t%0.12f\t%0.12f\n", (double)sum2 / counts2, (double)sum1 / counts1);
}

void write_cycle_rates(FILE *fp, fullerr_t *f)
{
	fputs("#Cycle\tRead 1 Error Rate\tRead 2 Error Rate\n", fp);
	for(uint64_t l = 0; l < f->l; ++l) {
		fprintf(fp, "%lu\t", l + 1);
		int sum1 = 0, sum2 = 0, counts1 = 0, counts2 = 0;
		for(int i = 0; i < 4; ++i) {
			sum1 += f->r1->qerr[i][l];
			counts1 += f->r1->qobs[i][l];
			sum2 += f->r2->qerr[i][l];
			counts2 += f->r2->qobs[i][l];
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
			if(j != nqscores - 1) {
				fprintf(ep, ","); fprintf(cp, ",");
            }
		}
		fprintf(ep, "|"); fprintf(cp, "|");
		for(j = 0; j < nqscores; ++j) {
			for(i = 0; i < 4; ++i) {
				fprintf(cp, i ? ":%lu": "%lu", f->r2->obs[i][j][l]);
				fprintf(ep, i ? ":%lu": "%lu", f->r2->err[i][j][l]);
			}
			if(j != nqscores - 1) {
				fprintf(ep, ","); fprintf(cp, ",");
            }
		}
		fprintf(ep, "\n"); fprintf(cp, "\n");
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

fullerr_t *fullerr_init(size_t l, char *bedpath, bam_hdr_t *hdr, int padding, int minFM, int maxFM, int flag) {
	fullerr_t *ret = (fullerr_t *)calloc(1, sizeof(fullerr_t));
	ret->l = l;
	ret->r1 = readerr_init(l);
	ret->r2 = readerr_init(l);
	if(bedpath) {
		ret->bed = kh_init(bed);
		ret->bed = parse_bed_hash(bedpath, hdr, padding);
	}
	ret->minFM = minFM;
	ret->maxFM = maxFM;
	ret->flag = flag;
	return ret;
}

void fullerr_destroy(fullerr_t *e) {
	if(e->r1) readerr_destroy(e->r1), e->r1 = NULL;
	if(e->r2) readerr_destroy(e->r2), e->r2 = NULL;
	cond_free(e->refcontig);
	if(e->bed) {
		kh_destroy(bed, e->bed);
	}
	free(e);
}

fmerr_t *fm_init(char *bedpath, bam_hdr_t *hdr, int padding, int flag) {
	fmerr_t *ret = (fmerr_t *)calloc(1, sizeof(fmerr_t *));
	if(bedpath && *bedpath) {
		ret->bed = kh_init(bed);
		ret->bed = parse_bed_hash(bedpath, hdr, padding);
	}
	ret->hash = kh_init(obs);
	ret->flag = flag;
	return ret;
}

void fm_destroy(fmerr_t *fm) {
	if(fm->bed) {
		kh_destroy(bed, fm->bed);
	}
	kh_destroy(obs, fm->hash);
	cond_free(fm->refcontig);
	free(fm);
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

	if(argc < 2) return err_usage(stderr, EXIT_FAILURE);

	if(strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0) err_usage(stderr, EXIT_SUCCESS);

	FILE *ofp = NULL, *d3 = NULL, *df = NULL, *dbc = NULL, *dc = NULL, *global_fp = NULL;
	char refcontig[200] = "";
	char *bedpath = NULL;
	int padding = -1;
	int minFM = 0;
	int maxFM = INT_MAX;
	int flag = 0;
	while ((c = getopt(argc, argv, "p:b:r:c:n:f:3:o:g:m:M:h?d")) >= 0) {
		switch (c) {
		case 'd': flag |= REQUIRE_DUPLEX; break;
		case 'm': minFM = atoi(optarg); break;
		case 'M': maxFM = atoi(optarg); break;
		case 'f': df = fopen(optarg, "w"); break;
		case 'o': strcpy(outpath, optarg); break;
		case '3': d3 = fopen(optarg, "w"); break;
		case 'c': dc = fopen(optarg, "w"); break;
		case 'n': dbc = fopen(optarg, "w"); break;
		case 'r': strcpy(refcontig, optarg); break;
		case 'b': bedpath = strdup(optarg); break;
		case 'p': padding = atoi(optarg); break;
		case 'g': global_fp = fopen(optarg, "w"); break;
		case '?': case 'h': return err_usage(stderr, EXIT_SUCCESS);
		}
	}

	if(padding < 0) {
		LOG_INFO("Padding not set. Setting to default value %i.\n", DEFAULT_PADDING);
	}

	if(!*outpath) {
		LOG_ERROR("Required -o parameter unset. Abort!\n");
	}
	ofp = open_ofp(outpath);

	if (argc != optind+2)
		return err_usage(stderr, EXIT_FAILURE);

	faidx_t *fai = fai_load(argv[optind]);

	fp = sam_open_format(argv[optind + 1], "r", &open_fmt);
	if (fp == NULL) {
		LOG_ERROR("Cannot open input file \"%s\"", argv[optind]);
	}

	header = sam_hdr_read(fp);
	if (header == NULL) {
		LOG_ERROR("Failed to read header for \"%s\"", argv[optind]);
	}
	// Get read length from the first read
	bam1_t *b = bam_init1();
	c = sam_read1(fp, header, b);
	fullerr_t *f = fullerr_init((size_t)b->core.l_qseq, bedpath, header, padding, minFM, maxFM, flag);
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
		write_3d_offsets(d3, f), fclose(d3), d3 = NULL;
	if(df)
		write_full_rates(df, f), fclose(df), df = NULL;
	if(dbc)
		write_base_rates(dbc, f), fclose(dbc), dbc = NULL;
	if(dc)
		write_cycle_rates(dc, f), fclose(dc), dc = NULL;
	if(!global_fp) global_fp = stdout;
	write_global_rates(global_fp, f);
	fclose(global_fp);
	fullerr_destroy(f);
	fclose(ofp);
	return EXIT_SUCCESS;
}
