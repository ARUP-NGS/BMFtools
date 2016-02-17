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
		if(!aux->iter) LOG_ERROR("Need to access bam with index.\n");
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

void process_pileup(bcf1_t *ret, const bam_pileup1_t *plp, int n_plp, int pos, int tid, vcfFile *outfp, bcf_hdr_t *hdr) {
	int khr, s, s2;
	khiter_t k;
	uint32_t *FA1, *PV1, *FA2, *PV2;
	std::string qname;
	uint8_t *seq, *seq2, *tmptag;
	// Build overlap hash
	std::unordered_map<std::string, BMF::UniqueObservation> obs;
	std::unordered_map<std::string, BMF::UniqueObservation>::const_iterator found;
	const int sk = 1;
	// Set the r1/r2 flags for the reads to ignore to 0
	// Set the ones where we see it twice to (BAM_FREAD1 | BAM_FREAD2).
	for(int i = 0; i < n_plp; ++i) {
		if(plp[i].is_del || plp[i].is_refskip) continue;
		// Skip any reads failed for FA < minFA or FR < minFR
		qname = std::string(bam_get_qname(plp[i].b));
		if((found = obs.find(qname)) == obs.end())
			obs[qname] = BMF::UniqueObservation(plp[i]);
		else ((BMF::UniqueObservation)(found->second)).add_obs(plp[i]);
	}
	BMF::VCFPos vcfline = BMF::VCFPos(obs, tid, pos, obs.size());
	vcfline.to_bcf(ret);
	bcf_write(outfp, hdr, ret);
}

int stack_main(int argc, char *argv[]) {
	return 0;
}
