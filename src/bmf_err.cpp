#include "bmf_err.h"


#define min_obs 10000uL

int err_main_main(int argc, char *argv[]);
int err_fm_main(int argc, char *argv[]);
int err_cycle_main(int argc, char *argv[]);

uint64_t get_max_obs(khash_t(obs) *hash)
{
	uint64_t ret = 0;
	for(khiter_t k = kh_begin(hash); k != kh_end(hash); ++k)
		if(kh_exist(hash, k))
			if(kh_val(hash, k).obs > ret)
				ret = kh_val(hash, k).obs;
	return ret;
}


uint64_t get_max_err(khash_t(obs) *hash)
{
	uint64_t ret = 0;
	for(khiter_t k = kh_begin(hash); k != kh_end(hash); ++k)
		if(kh_exist(hash, k))
			if(kh_val(hash, k).err > ret)
				ret = kh_val(hash, k).err;
	return ret;
}

int err_main_usage(FILE *fp, int retcode)
{
	fprintf(fp,
			"Usage: bmftools err main -o <out.tsv> <reference.fasta> <input.csrt.bam>\n"
			"Flags:\n"
			"-h/-?\t\tThis helpful help menu!\n"
			"-o\t\tPath to output file. Set to '-' or 'stdout' to emit to stdout.\n"
			"-a\t\tSet minimum mapping quality for inclusion.\n"
			"-$\t\tSet minimum calculated PV tag value for inclusion.\n"
			"-r:\t\tName of contig. If set, only reads aligned to this contig are considered\n"
			"-3:\t\tPath to write the 3d offset array in tabular format.\n"
			"-f:\t\tPath to write the full measured error rates in tabular format.\n"
			"-n:\t\tPath to write the cycle/nucleotide call error rates in tabular format.\n"
			"-c:\t\tPath to write the cycle error rates in tabular format.\n"
			"-b:\t\tPath to bed file for restricting analysis.\n"
			"-m:\t\tMinimum family size for inclusion. Default: 0.\n"
			"-M:\t\tMaximum family size for inclusion. Default: %i.\n"
			"-d:\t\tFlag to only calculate error rates for duplex reads.\n"
			"-D:\t\tFlag to only calculate error rates for non-duplex reads.\n"
			"-p:\t\tSet padding for bed region. Default: 50.\n"
			"-P:\t\tOnly include proper pairs.\n"
			, INT_MAX);
	exit(retcode);
	return retcode;
}

int err_fm_usage(FILE *fp, int retcode)
{
	fprintf(fp,
			"Usage: bmftools err fm -o <out.tsv> <reference.fasta> <input.csrt.bam>\n"
			"Flags:\n"
			"-h/-?\t\tThis helpful help menu!\n"
			"-o\t\tPath to output file. Set to '-' or 'stdout' to emit to stdout.\n"
			"-$\t\tSet minimum calculated PV tag value for inclusion.\n"
			"-a\t\tSet minimum mapping quality for inclusion.\n"
			"-r:\t\tName of contig. If set, only reads aligned to this contig are considered\n"
			"-b:\t\tPath to bed file for restricting analysis.\n"
			"-d:\t\tFlag to only calculate error rates for duplex reads.\n"
			"-p:\t\tSet padding for bed region. Default: 50.\n"
			"-P:\t\tOnly include proper pairs.\n"
			"-F:\t\tRequire that the FP tag be present and nonzero.\n"
			);
	exit(retcode);
	return retcode; // This never happens.
}


int err_cycle_usage(FILE *fp, int retcode)
{
	fprintf(fp,
			"Usage: bmftools err cycle -o <out.tsv> <reference.fasta> <input.csrt.bam>\n"
			"Flags:\n"
			"-h/-?\t\tThis helpful help menu!\n"
			"-o\t\tPath to output file. Set to '-' or 'stdout' to emit to stdout.\n"
			"-a\t\tSet minimum mapping quality for inclusion.\n"
			"-r:\t\tName of contig. If set, only reads aligned to this contig are considered\n"
			"-b:\t\tPath to bed file for restricting analysis.\n"
			"-p:\t\tSet padding for bed region. Default: 50.\n"
			"-P:\t\tOnly include proper pairs.\n"
			);
	exit(retcode);
	return retcode; // This never happens.
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

void err_cycle_report(FILE *fp, cycle_err_t *ce)
{
	LOG_DEBUG("Beginning err cycle report.\n");
	fprintf(fp, "#Cycle number\tRead 1 error\tRead 2 error\tRead 1 error counts"
			"\tRead 1 observation counts\tRead 2 error counts\tRead 2 observation counts\n");
	for(int i = 0; i < ce->rlen; ++i) {
		fprintf(fp, "%i\t%0.12f\t%lu\t%lu\t%0.12f\t%lu\t%lu\n",
				i + 1,
				(double)ce->r1[i].err / ce->r1[i].obs,
				ce->r1[i].err, ce->r1[i].obs,
				(double)ce->r2[i].err / ce->r2[i].obs, ce->r2[i].err, ce->r2[i].obs);
	}
}

int int32_cmp(const void *a, const void *b)
{
	return *((int *)a) - *((int *)b);
}

void err_fm_report(FILE *fp, fmerr_t *f)
{
	LOG_DEBUG("Beginning err fm report.\n");
	int khr, fm;
	khiter_t k, k1, k2;
	// Make a set of all FMs to print out.
	khash_t(obs_union) *key_union = kh_init(obs_union);
	for(k1 = kh_begin(f->hash1); k1 != kh_end(f->hash1); ++k1)
		if(kh_exist(f->hash1, k1))
			if((k = kh_get(obs_union, key_union, kh_key(f->hash1, k1))) == kh_end(key_union))
				k = kh_put(obs_union, key_union, kh_key(f->hash1, k1), &khr);
	for(k2 = kh_begin(f->hash2); k2 != kh_end(f->hash2); ++k2)
		if(kh_exist(f->hash2, k2))
			if((k = kh_get(obs_union, key_union, kh_key(f->hash2, k2))) == kh_end(key_union))
				k = kh_put(obs_union, key_union, kh_key(f->hash2, k2), &khr);
#if !NDEBUG
	for(k1 = kh_begin(f->hash1); k1 != kh_end(f->hash1); ++k1)
		if(kh_exist(f->hash1, k1))
			fprintf(stderr, "Key for FM: %i.\n",  kh_key(f->hash1, k1));
#endif

	// Write  header
	fprintf(fp, "##PARAMETERS\n##refcontig:\"%s\"\n##bed:\"%s\"\n"
			"##minMQ:%i\n##Duplex Required: %s.\n##Duplex Refused: %s.\n", f->refcontig ? f->refcontig: "N/A",
			f->bedpath? f->bedpath: "N/A", f->minMQ,
			f->flag & REQUIRE_DUPLEX ? "True": "False",
			f->flag & REFUSE_DUPLEX ? "True": "False");
	fprintf(fp, "##STATS\n##nread:%lu\n##nskipped:%lu\n", f->nread, f->nskipped);
	fprintf(fp, "#FM\tRead 1 Error\tRead 2 Error\tRead 1 Errors\tRead 1 Counts\tRead 2 Errors\tRead 2 Counts\n");
	int *tmp = (int *)malloc(key_union->n_occupied * sizeof(int));
	for(k = kh_begin(key_union), khr = 0; k != kh_end(key_union); ++k)
		if(kh_exist(key_union, k))
			tmp[khr++] = kh_key(key_union, k);
	qsort(tmp, key_union->n_occupied, sizeof(int), &int32_cmp);
	for(unsigned i = 0; i < key_union->n_occupied; ++i) {
		fm = tmp[i];
		fprintf(fp, "%i\t", fm);

		if((k1 = kh_get(obs, f->hash1, fm)) == kh_end(f->hash1))
			fprintf(fp, "-nan\t");
		else
			fprintf(fp, "%0.12f\t", (double)kh_val(f->hash1, k1).err / kh_val(f->hash1, k1).obs);
		LOG_DEBUG("R1 FM %i err, obs: %lu, %lu.\n", kh_key(f->hash1, k1), kh_val(f->hash1, k1).err, kh_val(f->hash1, k1).obs);

		if((k2 = kh_get(obs, f->hash2, fm)) == kh_end(f->hash2))
			fprintf(fp, "-nan\t");
		else
			fprintf(fp, "%0.12f\t", (double)kh_val(f->hash2, k2).err / kh_val(f->hash2, k2).obs);

		if(k1 != kh_end(f->hash1))
			fprintf(fp, "%lu\t%lu\t",
				kh_val(f->hash1, k1).err, kh_val(f->hash1, k1).obs);
		else fputs("0\t0\t", fp);
		if(k2 != kh_end(f->hash2))
			fprintf(fp, "%lu\t%lu\n",
				kh_val(f->hash2, k2).err, kh_val(f->hash2, k2).obs);
		else fputs("0\t0\n", fp);
	}
	free(tmp);
	kh_destroy(obs_union, key_union);
}

void err_report(FILE *fp, fullerr_t *e)
{
	LOG_DEBUG("Beginning error main report.\n");
	fprintf(fp, "{\n{\"total_read\": %lu},\n{\"total_skipped\": %lu},\n", e->nread, e->nskipped);
	uint64_t n1_obs = 0, n1_err = 0, n1_ins = 0;
	uint64_t n2_obs = 0, n2_err = 0, n2_ins = 0;
	// n_ins is number with insufficient observations to report.
	for(int i = 0; i < 4; ++i) {
		for(unsigned j = 0; j < nqscores; ++j) {
			for(unsigned k = 0; k < e->l; ++k) {
				n1_obs += e->r1->obs[i][j][k]; n1_err += e->r1->err[i][j][k];
				n2_obs += e->r2->obs[i][j][k]; n2_err += e->r2->err[i][j][k];
				if(e->r1->obs[i][j][k] < min_obs) ++n1_ins;
				if(e->r2->obs[i][j][k] < min_obs) ++n2_ins;
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
		for(unsigned j = 0; j < nqscores; ++j) {
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


void err_fm_core(char *fname, faidx_t *fai, fmerr_t *f, htsFormat *open_fmt)
{
	samFile *fp = sam_open_format(fname, "r", open_fmt);
	bam_hdr_t *hdr = sam_hdr_read(fp);
	if (!hdr) {
		LOG_ERROR("Failed to read input header from bam %s. Abort!\n", fname);
	}
	int32_t is_rev, ind, s, i, fc, rc, r, khr, DR, FP, FM, reflen, length, pos, tid_to_study = -1, last_tid = -1;
	char *ref = NULL; // Will hold the sequence for a  chromosome
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
	khash_t(obs) *hash;
	uint8_t *seq, *drdata, *fpdata;
	uint32_t *cigar, *pv_array;
	khiter_t k;
	bam1_t *b = bam_init1();
	while(LIKELY((r = sam_read1(fp, hdr, b)) != -1)) {
		if(++f->nread % 1000000 == 0) {
			LOG_INFO("Records read: %lu.\n", f->nread);
		}
		pv_array = (uint32_t *)array_tag(b, "PV");
		FM = bam_aux2i(bam_aux_get(b, "FM"));
		DR = ((drdata = bam_aux_get(b, "DR")) != NULL) ? bam_aux2i(drdata): 0;
		FP = ((fpdata = bam_aux_get(b, "FP")) != NULL) ? bam_aux2i(fpdata): 0;
		// Pass reads without FP tag.
		if((b->core.flag & 772) || // UNMAPPED, SECONDARY, QCFAIL
				b->core.qual < f->minMQ || //  minMQ
				(f->refcontig && tid_to_study != b->core.tid) ||
				((f->flag & REQUIRE_PROPER) && (!(b->core.flag & BAM_FPROPER_PAIR))) || // improper pair
			(f->bed && bed_test(b, f->bed) == 0) || // Outside of  region
			((f->flag & REQUIRE_DUPLEX) && !DR) || // Requires  duplex
			((f->flag & REFUSE_DUPLEX) && DR) ||
			((f->flag & REQUIRE_FP_PASS) && FP == 0)) /* Fails barcode QC */ {
				++f->nskipped;
				//LOG_DEBUG("Skipped record with name %s.\n", bam_get_qname(b));
				continue;
		}
		seq = (uint8_t *)bam_get_seq(b);
		cigar = bam_get_cigar(b);
#if !NDEBUG
		ifn_abort(cigar);
		ifn_abort(seq);
#endif
		hash = (b->core.flag & BAM_FREAD1) ? f->hash1: f->hash2;
		if(b->core.tid != last_tid) {
			cond_free(ref);
			LOG_DEBUG("Loading ref sequence for contig with name %s.\n", hdr->target_name[b->core.tid]);
			ref = fai_fetch(fai, hdr->target_name[b->core.tid], &reflen);
			if(!ref) {
				LOG_ERROR("Failed to load ref sequence for contig '%s'. Abort!\n", hdr->target_name[b->core.tid]);
			}
			last_tid = b->core.tid;
		}
		pos = b->core.pos;
		is_rev = b->core.flag & BAM_FREVERSE;
		//LOG_DEBUG("Max err count in hashmap: %lu.\n", get_max_err(hash));
		//LOG_DEBUG("Max obs count in hashmap: %lu.\n", get_max_obs(hash));
		if((k = kh_get(obs, hash, FM)) == kh_end(hash)) {
			k = kh_put(obs, hash, FM, &khr);
			memset(&kh_val(hash, k), 0, sizeof(obserr_t));
		}
		for(i = 0, rc = 0, fc = 0; i < b->core.n_cigar; ++i) {
			length = bam_cigar_oplen(*cigar);
			switch(bam_cigar_op(*cigar++)) {
			case BAM_CMATCH:
			case BAM_CEQUAL:
			case BAM_CDIFF:
				for(ind = 0; ind < length; ++ind) {
						if(pv_array[is_rev ? b->core.l_qseq - 1 - ind - rc: ind + rc] < f->minPV)
							continue;
					s = bam_seqi(seq, ind + rc);
					//fprintf(stderr, "Bi value: %i. s: %i.\n", bi, s);
					if(s == HTS_N || ref[pos + fc + ind] == 'N') continue;
					++kh_val(hash, k).obs;
					if(seq_nt16_table[(int8_t)ref[pos + fc + ind]] != s)
						++kh_val(hash, k).err;
				}
				rc += length; fc += length;
				break;
			case BAM_CSOFT_CLIP:
			case BAM_CHARD_CLIP:
			case BAM_CINS:
				rc += length;
				break;
			case BAM_CREF_SKIP:
			case BAM_CDEL:
				fc += length;
				break;
			}
		}
	}
	LOG_INFO("Total records read: %lu. Total records skipped: %lu.\n", f->nread, f->nskipped);
	cond_free(ref);
	bam_destroy1(b);
	bam_hdr_destroy(hdr), sam_close(fp);
}



cycle_err_t *err_cycle_core(char *fname, faidx_t *fai, htsFormat *open_fmt,
							char *bedpath, char *refcontig, int padding, unsigned minMQ, int flag)
{
	bam1_t *b = bam_init1();
	samFile *fp = sam_open_format(fname, "r", open_fmt);
	bam_hdr_t *hdr = sam_hdr_read(fp);
	if (!hdr) {
		LOG_ERROR("Failed to read input header from bam %s. Abort!\n", fname);
	}
	if(sam_read1(fp, hdr, b) == -1) {
		LOG_ERROR("Could not read bam record from bam %s. Abort!\n", fname);
	}
	const int32_t rlen = b->core.l_qseq;
	cycle_err_t *ce = cycle_init(bedpath, hdr, refcontig, padding, minMQ, rlen, flag);
	int32_t is_rev, ind, s, i, fc, rc, r, reflen, length, cycle, pos, tid_to_study = -1, last_tid = -1;
	char *ref = NULL; // Will hold the sequence for a  chromosome
	if(ce->refcontig) {
		for(int i = 0; i < hdr->n_targets; ++i) {
			if(!strcmp(hdr->target_name[i], ce->refcontig)) {
				tid_to_study = i; break;
			}
		}
		if(tid_to_study < 0) {
			LOG_ERROR("Contig %s not found in bam header. Abort mission!\n", ce->refcontig);
		}
	}
	uint8_t *seq, *fpdata;
	uint32_t *cigar;
	obserr_t *arr;
	cycle_loop(ce, b, seq, cigar, last_tid, ref, hdr, pos, is_rev, length, i, rc, fc, ind, s, cycle, arr, rlen, fpdata);
	while(LIKELY((r = sam_read1(fp, hdr, b)) != -1)) {
		cycle_loop(ce, b, seq, cigar, last_tid, ref, hdr, pos, is_rev, length, i, rc, fc, ind, s, cycle, arr, rlen, fpdata);
	}
	LOG_INFO("Total records read: %lu. Total records skipped: %lu.\n", ce->nread, ce->nskipped);
	cond_free(ref);
	bam_destroy1(b);
	bam_hdr_destroy(hdr), sam_close(fp);
	return ce;
}

void err_main_core(char *fname, faidx_t *fai, fullerr_t *f, htsFormat *open_fmt)
{
	if(!f->r1) f->r1 = readerr_init(f->l);
	if(!f->r2) f->r2 = readerr_init(f->l);
	samFile *fp = sam_open_format(fname, "r", open_fmt);
	bam_hdr_t *hdr = sam_hdr_read(fp);
	if (!hdr) {
		LOG_ERROR("Failed to read input header from bam %s. Abort!\n", fname);
	}
	int32_t i, s, c, len, pos, FM, RV, rc, fc, last_tid = -1, tid_to_study = -1, cycle, is_rev;
	unsigned ind;
	bam1_t *b = bam_init1();
	char *ref = NULL; // Will hold the sequence for a  chromosome
	if(f->refcontig) {
		for(i = 0; i < hdr->n_targets; ++i) {
			if(!strcmp(hdr->target_name[i], f->refcontig)) {
				tid_to_study = i; break;
			}
		}
		if(tid_to_study < 0) {
			LOG_ERROR("Contig %s not found in bam header. Abort mission!\n", f->refcontig);
		}
	}
	uint8_t *fdata, *rdata, *pdata, *seq, *qual;
	uint32_t *cigar, *pv_array, length;
	readerr_t *r;
	while(LIKELY((c = sam_read1(fp, hdr, b)) != -1)) {
		fdata = bam_aux_get(b, "FM");
		rdata = bam_aux_get(b, "RV");
		pdata = bam_aux_get(b, "FP");
		FM = fdata ? bam_aux2i(fdata): 0;
		RV = rdata ? bam_aux2i(rdata): 0;
		if((b->core.flag & 772) || b->core.qual < f->minMQ || (f->refcontig && tid_to_study != b->core.tid) ||
			(f->bed && bed_test(b, f->bed) == 0) || // Outside of region
			(FM < f->minFM) || (FM > f->maxFM) || // minFM 
			((f->flag & REQUIRE_PROPER) && (!(b->core.flag & BAM_FPROPER_PAIR))) || // skip improper pairs
			((f->flag & REQUIRE_DUPLEX) ? (RV == FM || RV == 0): ((f->flag & REFUSE_DUPLEX) && (RV != FM && RV != 0))) || // Requires 
			((f->flag & REQUIRE_FP_PASS) && pdata && bam_aux2i(pdata) == 0) /* Fails barcode QC */) {
				++f->nskipped;
				LOG_DEBUG("Skipped record with name %s.\n", bam_get_qname(b));
				continue;
		} // UNMAPPED, SECONDARY, QCFAIL, [removed SUPPLEMENTARY flag for a second],
		seq = (uint8_t *)bam_get_seq(b);
		qual = (uint8_t *)bam_get_qual(b);
		cigar = bam_get_cigar(b);
#if !NDEBUG
		ifn_abort(cigar);
		ifn_abort(seq);
		ifn_abort(qual);
#endif

		if(++f->nread % 1000000 == 0) {
			LOG_INFO("Records read: %lu.\n", f->nread);
		}
		if(b->core.tid != last_tid) {
			cond_free(ref);
			LOG_DEBUG("Loading ref sequence for contig with name %s.\n", hdr->target_name[b->core.tid]);
			ref = fai_fetch(fai, hdr->target_name[b->core.tid], &len);
			if(!ref) {
				LOG_ERROR("[Failed to load ref sequence for contig '%s'. Abort!\n", hdr->target_name[b->core.tid]);
			}
			last_tid = b->core.tid;
		}
		r = (b->core.flag & BAM_FREAD1) ? f->r1: f->r2;
		pos = b->core.pos;
		is_rev = (b->core.flag & BAM_FREVERSE);
		if(f->minPV) {
			pv_array = (uint32_t *)array_tag(b, "PV");
			for(i = 0, rc = 0, fc = 0; i < b->core.n_cigar; ++i) {
				length = bam_cigar_oplen(*cigar);
				switch(bam_cigar_op(*cigar++)) {
				case BAM_CMATCH:
				case BAM_CEQUAL:
				case BAM_CDIFF:
					for(ind = 0; ind < length; ++ind) {
						if(pv_array[is_rev ? b->core.l_qseq - 1 - ind - rc: ind + rc] < f->minPV)
							continue;
						s = bam_seqi(seq, ind + rc);
						//fprintf(stderr, "Bi value: %i. s: %i.\n", bi, s);
						if(s == HTS_N || ref[pos + fc + ind] == 'N') continue;
#if !NDEBUG
						if(UNLIKELY(qual[ind + rc] > nqscores + 1)) { // nqscores + 2 - 1
							LOG_ERROR("Quality score is too high. int: %i. char: %c. Max permitted: %lu.\n",
									(int)qual[ind + rc], qual[ind + rc], nqscores + 1);
						}
#endif
						cycle = is_rev ? b->core.l_qseq - 1 - ind - rc: ind + rc;
						++r->obs[bamseq2i[s]][qual[ind + rc] - 2][cycle];
						if(seq_nt16_table[(int8_t)ref[pos + fc + ind]] != s)
							++r->err[bamseq2i[s]][qual[ind + rc] - 2][cycle];
					}
					rc += length; fc += length;
					break;
				case BAM_CSOFT_CLIP:
				case BAM_CHARD_CLIP:
				case BAM_CINS:
					rc += length;
					break;
				case BAM_CREF_SKIP:
				case BAM_CDEL:
					fc += length;
					break;
				}
			}
		} else {
			for(i = 0, rc = 0, fc = 0; i < b->core.n_cigar; ++i) {
				length = bam_cigar_oplen(*cigar);
				switch(bam_cigar_op(*cigar++)) {
				case BAM_CMATCH:
				case BAM_CEQUAL:
				case BAM_CDIFF:
					for(ind = 0; ind < length; ++ind) {
						//LOG_DEBUG("Base at rpos %i, fpos %i: %c. Ref base: %c.\n", ind + rc + 1, pos + fc + ind + 1, seq_nt16_str[bam_seqi(seq, ind + rc)], ref[pos + fc + ind]);
						s = bam_seqi(seq, ind + rc);
						//fprintf(stderr, "Bi value: %i. s: %i.\n", bi, s);
						if(s == HTS_N || ref[pos + fc + ind] == 'N') continue;
#if !NDEBUG
						if(UNLIKELY(qual[ind + rc] > nqscores + 1)) { // nqscores + 2 - 1
							LOG_ERROR("Quality score is too high. int: %i. char: %c. Max permitted: %lu.\n",
									(int)qual[ind + rc], qual[ind + rc], nqscores + 1);
						}
#endif
						cycle = is_rev ? b->core.l_qseq - 1 - ind - rc: ind + rc;
						++r->obs[bamseq2i[s]][qual[ind + rc] - 2][cycle];
						if(seq_nt16_table[(int8_t)ref[pos + fc + ind]] != s) {
							//LOG_DEBUG("Hey, error found at rpos %i.\n", ind + rc + 1);
							++r->err[bamseq2i[s]][qual[ind + rc] - 2][cycle];
							//LOG_DEBUG("Error total at cycle %i : %lu.\n", cycle + 1, r->err[bamseq2i[s]][qual[ind + rc] - 2][cycle]);
						}
					}
					rc += length; fc += length;
					break;
				case BAM_CSOFT_CLIP:
				case BAM_CHARD_CLIP:
				case BAM_CINS:
					rc += length;
					break;
				case BAM_CREF_SKIP:
				case BAM_CDEL:
					fc += length;
					break;
				}
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
	char *ref = NULL; // Will hold the sequence for a contig
	int tid_to_study = -1;
	unsigned ind;
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
		if((b->core.flag & 772) || (f->refcontig && tid_to_study != b->core.tid)) {++f->nskipped; continue;} // UNMAPPED, SECONDARY, SUPPLEMENTARY, QCFAIL
		const uint8_t *seq = (uint8_t *)bam_get_seq(b);
		const uint8_t *qual = (uint8_t *)bam_get_qual(b);
		const uint32_t *cigar = bam_get_cigar(b);
#if !NDEBUG
		ifn_abort(cigar);
		ifn_abort(seq);
		ifn_abort(qual);
#endif

		if(UNLIKELY(++f->nread % 1000000 == 0)) {
			LOG_INFO("Records read: %lu.\n", f->nread);
		}
		if(b->core.tid != last_tid) {
			cond_free(ref);
			LOG_INFO("Loading ref sequence for contig with name %s.\n", hdr->target_name[b->core.tid]);
			ref = fai_fetch(fai, hdr->target_name[b->core.tid], &len);
			last_tid = b->core.tid;
		}
		// rc -> readcount
		// fc -> reference base  count
		//fprintf(stderr, "Pointer to readerr_t r: %p.\n", r);
		const int32_t pos = b->core.pos;
		for(int i = 0, rc = 0, fc = 0; i < b->core.n_cigar; ++i) {
			//fprintf(stderr, "Qual %p, seq %p, cigar %p.\n", seq, qual, cigar);
			int s; // seq value, base 
			const uint32_t len = bam_cigar_oplen(*cigar);
			switch(bam_cigar_op(*cigar++)) {
			case BAM_CMATCH:
			case BAM_CEQUAL:
			case BAM_CDIFF:
				for(ind = 0; ind < len; ++ind) {
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
	unsigned i, j;
	for(l = 0; l < f->l; ++l) {
		for(j = 0; j < nqscores; ++j) {
			for(i = 0; i < 4u; ++i) {
				if(f->r1->obs[i][j][l])
					fprintf(fp, i ? ":%0.12f": "%0.12f", (double)f->r1->err[i][j][l] / f->r1->obs[i][j][l]);
				else fputs(i ? ":-1337": "-1337", fp);
			}
			if(j != nqscores - 1) fputc(',', fp);
		}
		fputc('|', fp);
		for(j = 0; j < nqscores; ++j) {
			for(i = 0; i < 4u; ++i) {
				if(f->r2->obs[i][j][l])
					fprintf(fp, i ? ":%0.12f": "%0.12f", (double)f->r2->err[i][j][l] / f->r2->obs[i][j][l]);
				else fputs(i ? ":-1337": "-1337", fp);
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
			LOG_DEBUG("obs: %lu. err: %lu.\n", f->r1->qerr[i][l], f->r1->qobs[i][l]);
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
	uint64_t sum1 = 0, sum2 = 0, counts1 = 0, counts2 = 0;
	for(uint64_t l = 0; l < f->l; ++l) {
		for(int i = 0; i < 4; ++i) {
			sum1 += f->r1->qerr[i][l];
			counts1 += f->r1->qobs[i][l];
			sum2 += f->r2->qerr[i][l];
			counts2 += f->r2->qobs[i][l];
		}
	}
	fprintf(fp, "#Global Error Rates\t%0.12f\t%0.12f\n", (double)sum1 / counts1, (double)sum2 / counts2);
	fprintf(fp, "#Global Sum/Count1 Sum/Count2\t%lu\t%lu\t%lu\t%lu\n", sum1, counts1, sum2, counts2);
}

void write_cycle_rates(FILE *fp, fullerr_t *f)
{
	fputs("#Cycle\tRead 1 Error Rate\tRead 2 Error Rate\tRead 1 Error Count\t"
			"Read 1 Obs Count\tRead 2 Error Count\tRead 2 Obs Count\n", fp);
	for(uint64_t l = 0; l < f->l; ++l) {
		fprintf(fp, "%lu\t", l + 1);
		uint64_t sum1 = 0, sum2 = 0, counts1 = 0, counts2 = 0;
		for(int i = 0; i < 4; ++i) {
			sum1 += f->r1->qerr[i][l];
			counts1 += f->r1->qobs[i][l];
			sum2 += f->r2->qerr[i][l];
			counts2 += f->r2->qobs[i][l];
		}
		fprintf(fp, "%0.12f\t%0.12f\t%lu\t%lu\t%lu\t%lu\n", (double)sum1 / counts1, (double)sum2 / counts2,
				sum1, counts1, sum2, counts2);
	}
}

void impute_scores(fullerr_t *f)
{
	unsigned j, i;
	uint64_t l;
	for(i = 0; i < 4u; ++i)
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
			for(unsigned j = 1; j < nqscores; ++j) { // Skip qualities of 2
				f->r1->qpvsum[i][l] += pow(10., (double)(-0.1 * (j + 2))) * f->r1->obs[i][j][l];
				f->r2->qpvsum[i][l] += pow(10., (double)(-0.1 * (j + 2))) * f->r2->obs[i][j][l];
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
	for(int i = 0; i < 4; ++i) {
		for(unsigned j = 0; j < nqscores; ++j) {
			for(uint64_t l = 0; l < f->l; ++l) {
				if(f->r1->obs[i][j][l] > min_obs)
					f->r1->final[i][j][l] = pv2ph((double)f->r1->err[i][j][l] / f->r1->obs[i][j][l]);
				if(f->r2->obs[i][j][l] > min_obs)
					f->r2->final[i][j][l] = pv2ph((double)f->r2->err[i][j][l] / f->r2->obs[i][j][l]);
			}
		}
	}
}

void write_counts(fullerr_t *f, FILE *cp, FILE *ep)
{
	FILE *dictwrite = fopen("dict.txt", "w");
	fprintf(dictwrite, "{\n\t");
	unsigned i, j, l;
	for(l = 0; l < f->l; ++l) {
		for(j = 0; j < nqscores; ++j) {
			for(i = 0; i < 4u; ++i) {
				fprintf(dictwrite, "'r1,%c,%i,%u,obs': %lu,\n\t", NUM2NUC_STR[i], j + 2, l + 1, f->r1->obs[i][j][l]);
				fprintf(dictwrite, "'r2,%c,%i,%u,obs': %lu,\n\t", NUM2NUC_STR[i], j + 2, l + 1, f->r2->obs[i][j][l]);
				fprintf(dictwrite, "'r1,%c,%i,%u,err': %lu,\n\t", NUM2NUC_STR[i], j + 2, l + 1, f->r1->err[i][j][l]);
				if(i == 3 && j == nqscores - 1 && l == f->l - 1)
					fprintf(dictwrite, "'r2,%c,%i,%u,err': %lu\n}", NUM2NUC_STR[i], j + 2, l + 1, f->r2->err[i][j][l]);
				else
					fprintf(dictwrite, "'r2,%c,%i,%u,err': %lu,\n\t", NUM2NUC_STR[i], j + 2, l + 1, f->r2->err[i][j][l]);
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

fullerr_t *fullerr_init(size_t l, char *bedpath, bam_hdr_t *hdr, int padding, int minFM, int maxFM, int flag, int minMQ, uint32_t minPV) {
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
	ret->minMQ = minMQ;
	ret->minPV = minPV;
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

fmerr_t *fm_init(char *bedpath, bam_hdr_t *hdr, char *refcontig, int padding, int flag, int minMQ, uint32_t minPV) {
	fmerr_t *ret = (fmerr_t *)calloc(1, sizeof(fmerr_t));
	if(bedpath && *bedpath) {
		ret->bed = parse_bed_hash(bedpath, hdr, padding);
		ret->bedpath = strdup(bedpath);
	}
	if(refcontig && *refcontig) {
		ret->refcontig = strdup(refcontig);
	}
	ret->hash1 = kh_init(obs);
	ret->hash2 = kh_init(obs);
	ret->flag = flag;
	ret->minMQ = minMQ;
	return ret;
}

void fm_destroy(fmerr_t *fm) {
	if(fm->bed) kh_destroy(bed, fm->bed);
	kh_destroy(obs, fm->hash1);
	kh_destroy(obs, fm->hash2);
	cond_free(fm->refcontig);
	cond_free(fm->bedpath);
	free(fm);
}
/*
typedef struct cycle_err {
	int32_t flag; //
	uint32_t rlen;
	obserr_t *r1;
	obserr_t *r2;
	char *refcontig;
	char *bedpath;
} cycle_err_t; */

cycle_err_t *cycle_init(char *bedpath, bam_hdr_t *hdr, char *refcontig, int padding, int minMQ, int rlen, int flag)
{
	cycle_err_t *ret = (cycle_err_t *)calloc(1, sizeof(cycle_err_t));
	if(bedpath && *bedpath) {
		ret->bed = parse_bed_hash(bedpath, hdr, padding);
		ret->bedpath = strdup(bedpath);
	}
	if(refcontig && *refcontig) {
		ret->refcontig = strdup(refcontig);
	}
	ret->rlen = rlen;
	ret->r1 = (obserr_t *)calloc(rlen, sizeof(obserr_t));
	ret->r2 = (obserr_t *)calloc(rlen, sizeof(obserr_t));
	ret->minMQ = minMQ;
	ret->flag = flag;
	return ret;
}

void cycle_destroy(cycle_err_t *c)
{
	if(c->bed) kh_destroy(bed, c->bed), c->bed = NULL;
	cond_free(c->bedpath);
	cond_free(c->refcontig);
	cond_free(c->r1); cond_free(c->r2);
	free(c);
}


int err_usage(FILE *fp, int retcode)
{
	fprintf(stderr, "bmftools err\nSubcommands: main, fm, cycle.\n");
	exit(retcode);
	return retcode; // This neverhappens
}

int err_main(int argc, char *argv[])
{
	if(argc < 2) return err_usage(stderr, EXIT_FAILURE);
	if(strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0)
		return err_usage(stderr, EXIT_SUCCESS);
	if(strcmp(argv[1], "main") == 0)
		return err_main_main(argc - 1, argv + 1);
	if(strcmp(argv[1], "fm") == 0)
		return err_fm_main(argc - 1, argv + 1);
	if(strcmp(argv[1], "cycle") == 0)
		return err_cycle_main(argc - 1, argv + 1);
	LOG_ERROR("Unrecognized subcommand '%s'. Abort!\n", argv[1]);
}


int err_main_main(int argc, char *argv[])
{
	htsFormat open_fmt = (htsFormat){sequence_data, bam, {1, 3}, gzip, 0, NULL};
	samFile *fp = NULL;
	bam_hdr_t *header = NULL;
	int c, minMQ = 0;
	char outpath[500] = "";

	if(argc < 2) return err_main_usage(stderr, EXIT_FAILURE);

	if(strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0) err_main_usage(stderr, EXIT_SUCCESS);



	FILE *ofp = NULL, *d3 = NULL, *df = NULL, *dbc = NULL, *dc = NULL, *global_fp = NULL;
	char refcontig[200] = "";
	char *bedpath = NULL;
	int padding = -1;
	int minFM = 0;
	int maxFM = INT_MAX;
	int flag = 0;
	uint32_t minPV = 0;
	while ((c = getopt(argc, argv, "a:p:b:r:c:n:f:3:o:g:m:M:$:h?FdDP")) >= 0) {
		switch (c) {
		case 'a': minMQ = atoi(optarg); break;
		case 'd': flag |= REQUIRE_DUPLEX; break;
		case 'D': flag |= REFUSE_DUPLEX; break;
		case 'P': flag |= REQUIRE_PROPER; break;
		case 'F': flag |= REQUIRE_FP_PASS; break;
		case 'm': minFM = atoi(optarg); break;
		case 'M': maxFM = atoi(optarg); break;
		case 'f': df = open_ofp(optarg); break;
		case 'o': strcpy(outpath, optarg); break;
		case '3': d3 = open_ofp(optarg); break;
		case 'c': dc = open_ofp(optarg); break;
		case 'n': dbc = open_ofp(optarg); break;
		case 'r': strcpy(refcontig, optarg); break;
		case 'b': bedpath = strdup(optarg); break;
		case 'p': padding = atoi(optarg); break;
		case 'g': global_fp = open_ofp(optarg); break;
		case '$': minPV = strtoul(optarg, NULL, 0); break;
		case '?': case 'h': return err_main_usage(stderr, EXIT_SUCCESS);
		}
	}

	if(padding < 0 && bedpath && *bedpath) {
		LOG_INFO("Padding not set. Setting to default value %i.\n", DEFAULT_PADDING);
	}

	ofp = *outpath ? open_ofp(outpath): NULL;

	if (argc != optind+2)
		return err_main_usage(stderr, EXIT_FAILURE);

	faidx_t *fai = fai_load(argv[optind]);

	if ((fp = sam_open_format(argv[optind + 1], "r", &open_fmt)) == NULL) {
		LOG_ERROR("Cannot open input file \"%s\"", argv[optind]);
	}
	if ((header = sam_hdr_read(fp)) == NULL) {
		LOG_ERROR("Failed to read header for \"%s\"", argv[optind]);
	}

	if(minPV) check_bam_tag_exit(argv[optind + 1], "PV");
	if(minFM || maxFM != INT_MAX) check_bam_tag_exit(argv[optind + 1], "FM");

	// Get read length from the first read
	bam1_t *b = bam_init1();
	c = sam_read1(fp, header, b);
	fullerr_t *f = fullerr_init((size_t)b->core.l_qseq, bedpath, header,
								 padding, minFM, maxFM, flag, minMQ, minPV);
	sam_close(fp);
	fp = NULL;
	bam_destroy1(b);
	if(*refcontig) f->refcontig = strdup(refcontig);
	bam_hdr_destroy(header), header = NULL;
	err_main_core(argv[optind + 1], fai, f, &open_fmt);
	LOG_DEBUG("Core finished.\n");
	fai_destroy(fai);
	fill_qvals(f);
	impute_scores(f);
	fill_sufficient_obs(f);
	if(ofp)
		write_final(ofp, f), fclose(ofp);
	if(d3)
		write_3d_offsets(d3, f), fclose(d3), d3 = NULL;
	if(df)
		write_full_rates(df, f), fclose(df), df = NULL;
	if(dbc)
		write_base_rates(dbc, f), fclose(dbc), dbc = NULL;
	if(dc)
		write_cycle_rates(dc, f), fclose(dc), dc = NULL;
	if(!global_fp) global_fp = stdout;
	write_global_rates(global_fp, f); fclose(global_fp);
	fullerr_destroy(f);
	return EXIT_SUCCESS;
}

int err_cycle_main(int argc, char *argv[])
{
	htsFormat open_fmt = (htsFormat){sequence_data, bam, {1, 3}, gzip, 0, NULL};
	char outpath[500] = "";

	if(argc < 2) return err_cycle_usage(stderr, EXIT_FAILURE);

	if(strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0) return err_cycle_usage(stderr, EXIT_SUCCESS);



	FILE *ofp = NULL;
	char refcontig[200] = "";
	char *bedpath = NULL;
	int padding = -1, minMQ = 0, flag = 0, c;
	while ((c = getopt(argc, argv, "p:b:r:o:a:h?P")) >= 0) {
		switch (c) {
		case 'a': minMQ = atoi(optarg); break;
		case 'o': strcpy(outpath, optarg); break;
		case 'r': strcpy(refcontig, optarg); break;
		case 'b': bedpath = strdup(optarg); break;
		case 'p': padding = atoi(optarg); break;
		case 'P': flag |= REQUIRE_PROPER; break;
		case '?': case 'h': return err_cycle_usage(stderr, EXIT_SUCCESS);
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
		return err_cycle_usage(stderr, EXIT_FAILURE);

	faidx_t *fai = fai_load(argv[optind]);

	cycle_err_t *ce = err_cycle_core(argv[optind + 1], fai, &open_fmt, bedpath, refcontig, padding, minMQ, flag);
	err_cycle_report(ofp, ce); fclose(ofp);
	cycle_destroy(ce);
	fai_destroy(fai);
	return EXIT_SUCCESS;
}


int err_fm_main(int argc, char *argv[])
{
	htsFormat open_fmt;
	memset(&open_fmt, 0, sizeof(htsFormat));
	open_fmt.category = sequence_data;
	open_fmt.format = bam;
	open_fmt.version.major = 1;
	open_fmt.version.minor = 3;
	samFile *fp = NULL;
	bam_hdr_t *header = NULL;
	char outpath[500] = "";

	if(argc < 2) return err_fm_usage(stderr, EXIT_FAILURE);

	if(strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0) return err_fm_usage(stderr, EXIT_SUCCESS);



	FILE *ofp = NULL;
	char refcontig[200] = "";
	char *bedpath = NULL;
	int flag = 0, padding = -1, minMQ = 0, c;
	uint32_t minPV;
	while ((c = getopt(argc, argv, "$:p:b:r:o:a:Fh?dP")) >= 0) {
		switch (c) {
		case 'a': minMQ = atoi(optarg); break;
		case 'd': flag |= REQUIRE_DUPLEX; break;
		case 'o': strcpy(outpath, optarg); break;
		case 'r': strcpy(refcontig, optarg); break;
		case 'b': bedpath = strdup(optarg); break;
		case 'p': padding = atoi(optarg); break;
		case 'P': flag |= REQUIRE_PROPER; break;
		case '$': minPV = strtoul(optarg, NULL, 0); break;
		case 'F': flag |= REQUIRE_FP_PASS; break;
		case '?': case 'h': return err_fm_usage(stderr, EXIT_SUCCESS);
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
		return err_fm_usage(stderr, EXIT_FAILURE);

	faidx_t *fai = fai_load(argv[optind]);

	fp = sam_open_format(argv[optind + 1], "r", &open_fmt);
	if (fp == NULL) {
		LOG_ERROR("Cannot open input file \"%s\"", argv[optind]);
	}
	check_bam_tag_exit(argv[optind + 1], "FM");
	check_bam_tag_exit(argv[optind + 1], "FP");
	check_bam_tag_exit(argv[optind + 1], "RV");
	if(flag & (REQUIRE_DUPLEX | REFUSE_DUPLEX)) check_bam_tag_exit(argv[optind + 1], "DR");

	header = sam_hdr_read(fp);
	if (header == NULL) {
		LOG_ERROR("Failed to read header for \"%s\"", argv[optind]);
	}
	fmerr_t *f = fm_init(bedpath, header, refcontig, padding, flag, minMQ, minPV);
	// Get read length from the first 
	bam_hdr_destroy(header); header = NULL;
	err_fm_core(argv[optind + 1], fai, f, &open_fmt);
	err_fm_report(ofp, f); fclose(ofp);
	fai_destroy(fai);
	fm_destroy(f);
	return EXIT_SUCCESS;
}
