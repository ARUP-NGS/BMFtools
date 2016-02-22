#include "bmf_stack.h"
namespace {
	const char *vcf_header_lines[] =  {
			"##INFO=<ID=BMF_VET,Number=A,Type=Integer,Description=\"1 if the variant passes vetting, 0 otherwise.\">",
			"##INFO=<ID=BMF_UNIOBS,Number=A,Type=Integer,Description=\"Number of unique observations supporting variant at position.\">",
			"##INFO=<ID=BMF_DUPLEX,Number=A,Type=Integer,Description=\"Number of duplex reads supporting variant at position.\">",
			"##INFO=<ID=BMF_FAIL,Number=A,Type=Integer,Description=\"Number of unique observations at position failing filters.\">",
			"##INFO=<ID=DUPLEX_DEPTH,Number=1,Type=Integer,Description=\"Number of duplex reads passing filters at position.\">",
			"##INFO=<ID=DISC_OVERLAP,Number=1,Type=Integer,Description=\"Number of read pairs at position with discordant base calls.\">",
			"##INFO=<ID=OVERLAP,Number=1,Type=Integer,Description=\"Number of overlapping read pairs combined into single observations at position.\">"
	};
}

typedef struct {
	samFile *fp;
	hts_itr_t *iter;
	bam_hdr_t *header;
	vcfFile *vcf_fp;
	vcfFile *vcf_ofp;
	bcf_hdr_t *vcf_header;
	khash_t(bed) *bed;
	float minFR; // Minimum fraction of family members agreed on base
	float minAF; // Minimum aligned fraction
	int max_depth;
	int minFM;
	uint32_t minFA;
	uint32_t minPV;
	uint32_t minMQ;
	int minCount;
	int minDuplex;
	int minOverlap;
	int skip_improper;
	int vet_all; // If provided an unindexed variant file, vet all variants, not just those in bed region of interest.
	uint32_t skip_flag; // Skip reads with any bits set to true
} aux_t;

void stack_usage(int retcode)
{
	const char *buf =
			"Usage:\nbmftools vet -o <out.vcf [stdout]> <in.vcf.gz/in.bcf> <in.srt.indexed.bam>\n"
			"Optional arguments:\n"
			"-b, --bedpath\tPath to bed file to only validate variants in said region. REQUIRED.\n"
			"-c, --min-count\tMinimum number of observations for a given allele passing filters to pass variant. Default: 1.\n"
			"-s, --min-family-size\tMinimum number of reads in a family to include a that collapsed observation\n"
			"-2, --skip-secondary\tSkip secondary alignments.\n"
			"-S, --skip-supplementary\tSkip supplementary alignments.\n"
			"-q, --skip-qcfail\tSkip reads marked as QC fail.\n"
			"-f, --min-fraction-agreed\tMinimum fraction of reads in a family agreed on a base call\n"
			"-v, --min-phred-quality\tMinimum calculated p-value on a base call in phred space\n"
			"-p, --padding\tNumber of bases outside of bed region to pad.\n"
			"-a, --min-family-agreed\tMinimum number of reads in a family agreed on a base call\n"
			"-m, --min-mapping-quality\tMinimum mapping quality for reads for inclusion\n"
			"-B, --emit-bcf-format\tEmit bcf-formatted output. (Defaults to vcf).\n";
	fprintf(stderr, buf);
	exit(retcode);
}

static int read_bam(void *data, bam1_t *b)
{
	aux_t *aux = (aux_t*)data; // data in fact is a pointer to an auxiliary structure
	int ret;
	for(;;)
	{
		if(!aux->iter) LOG_EXIT("Need to access bam with index.\n");
		ret = sam_itr_next(aux->fp, aux->iter, b);
		if ( ret<0 ) break;
		// Skip unmapped, secondary, qcfail, duplicates.
		// Skip improper if option set
		// Skip MQ < minMQ
		// Skip FM < minFM
		// Skip AF < minAF
		if ((b->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP)) ||
			(aux->skip_improper && ((b->core.flag & BAM_FPROPER_PAIR) == 0)) || // Skip improper if set.
			(int)b->core.qual < aux->minMQ || (bam_aux2i(bam_aux_get(b, "FM")) < aux->minFM) ||
			(bam_aux2i(bam_aux_get(b, "FP")) == 0) || (aux->minAF && bam_aux2f(bam_aux_get(b, "AF")) < aux->minAF))
				continue;
		break;
	}
	return ret;
}

void process_pileup(const bam_pileup1_t *plp, int n_plp) {
	int khr, s, s2;
	khiter_t k;
	uint32_t *FA1, *PV1, *FA2, *PV2;
	char *qname;
	uint8_t *seq, *seq2, *tmptag, *drdata;
	// Build overlap hash
	khash_t(names) *hash = kh_init(names);
	const int sk = 1;
	// Set the r1/r2 flags for the reads to ignore to 0
	// Set the ones where we see it twice to (BAM_FREAD1 | BAM_FREAD2).
	for(int i = 0; i < n_plp; ++i) {
		if(plp[i].is_del || plp[i].is_refskip) continue;
		// Skip any reads failed for FA < minFA or FR < minFR
		qname = bam_get_qname(plp[i].b);
		k = kh_get(names, hash, qname);
		if(k == kh_end(hash)) {
			k = kh_put(names, hash, qname, &khr);
			kh_val(hash, k) = &plp[i];
		} else {
			bam_aux_append(plp[i].b, "SK", 'i', sizeof(int), (uint8_t *)&sk); // Skip
			bam_aux_append(kh_val(hash, k)->b, "KR", 'i', sizeof(int), (uint8_t *)&sk); // Keep Read
			if((tmptag = bam_aux_get(kh_val(hash, k)->b, "fm")) == NULL) {
				uint8_t *FM1 = bam_aux_get(kh_val(hash, k)->b, "FM");
				const int FM_sum = bam_aux2i(FM1) + bam_aux2i(bam_aux_get(plp[i].b, "FM"));
				bam_aux_del(kh_val(hash, k)->b, FM1);
				bam_aux_append(kh_val(hash, k)->b, "FM", 'i', sizeof(int), (uint8_t *)&FM_sum);
				bam_aux_append(kh_val(hash, k)->b, "fm", 'i', sizeof(int), (uint8_t *)&sk);
				bam_aux_append(plp[i].b, "fm", 'i', sizeof(int), (uint8_t *)&sk);
			}
			PV1 = (uint32_t *)array_tag(kh_val(hash, k)->b, "PV");
			FA1 = (uint32_t *)array_tag(kh_val(hash, k)->b, "FA");
			seq = bam_get_seq(kh_val(hash, k)->b);
			s = bam_seqi(seq, kh_val(hash, k)->qpos);
			PV2 = (uint32_t *)array_tag(plp[i].b, "PV");
			FA2 = (uint32_t *)array_tag(plp[i].b, "FA");
			seq2 = bam_get_seq(plp[i].b);
			s2 = bam_seqi(seq2, plp[i].qpos);
			const int32_t arr_qpos1 = arr_qpos(kh_val(hash, k));
			const int32_t arr_qpos2 = arr_qpos(&plp[i]);
			if(s == s2) {
				PV1[arr_qpos1] = agreed_pvalues(PV1[arr_qpos1], PV2[arr_qpos2]);
				FA1[arr_qpos1] = FA1[arr_qpos1] + FA2[arr_qpos2];
			} else if(s == HTS_N) {
				set_base(seq, seq_nt16_str[bam_seqi(seq2, plp[i].qpos)], kh_val(hash, k)->qpos);
				PV1[arr_qpos1] = PV2[arr_qpos2];
				FA1[arr_qpos1] = FA2[arr_qpos2];
			} else if(s2 != HTS_N) {
				// Disagreed, both aren't N: N the base, set agrees and p values to 0!
				n_base(seq, kh_val(hash, k)->qpos); // if s2 == HTS_N, do nothing.
				PV1[arr_qpos1] = 0u;
				FA1[arr_qpos1] = 0u;
			}
		}
	}
	for(int i = 0; i < n_plp; ++i) {
		if(plp[i].is_del || plp[i].is_refskip) continue;
		if((tmptag = bam_aux_get(plp[i].b, "SK")) != NULL) {
			continue;
		}

		seq = bam_get_seq(plp[i].b);
		FA1 = (uint32_t *)array_tag(plp[i].b, "FA");
		PV1 = (uint32_t *)array_tag(plp[i].b, "PV");
		bcf1_t *vrec = NULL;
		int j = 0;
		if(bam_seqi(seq, plp[i].qpos) == seq_nt16_table[(uint8_t)vrec->d.allele[j][0]]) { // Match!
			const int32_t arr_qpos1 = arr_qpos(&plp[i]);
			if(0) {
				continue;
			}
			if((drdata = bam_aux_get(plp[i].b, "DR")) != NULL && bam_aux2i(drdata)) {
			}
			if((tmptag = bam_aux_get(plp[i].b, "KR")) != NULL) {
				bam_aux_del(plp[i].b, tmptag);
			}
		}
	}
	for(int i = 0; i < n_plp; ++i)
		if((tmptag = bam_aux_get(plp[i].b, "SK")) != NULL) bam_aux_del(plp[i].b, tmptag);
	kh_destroy(names, hash);
}

int stack_main(int argc, char *argv[]) {
	return 0;
}
