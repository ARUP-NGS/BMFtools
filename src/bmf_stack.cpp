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




void stack_usage(int retcode)
{
	const char *buf =
			"Usage:\nbmftools stack -o <out.vcf [stdout]> <in.srt.indexed.bam>\n"
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


static int read_bam(aux_t *data, bam1_t *b)
{
	int ret;
	for(;;)
	{
		if(!data->iter) LOG_ERROR("Need to access bam with index.\n");
		ret = sam_itr_next(data->fp, data->iter, b);
		if ( ret<0 ) break;
		// Skip unmapped, secondary, qcfail, duplicates.
		// Skip improper if option set
		// Skip MQ < minMQ
		// Skip FM < minFM
		// Skip AF < minAF
		if ((b->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP)) ||
			(data->skip_improper && ((b->core.flag & BAM_FPROPER_PAIR) == 0)) || // Skip improper if set.
			(int)b->core.qual < data->minMQ || (bam_aux2i(bam_aux_get(b, "FM")) < data->minFM) ||
			(bam_aux2i(bam_aux_get(b, "FP")) == 0) || (data->minAF && bam_aux2f(bam_aux_get(b, "AF")) < data->minAF))
				continue;
		break;
	}
	return ret;
}

void process_pileup(bcf1_t *ret, const bam_pileup1_t *plp, int n_plp, int pos, int tid, vcfFile *outfp, bcf_hdr_t *hdr) {
	std::string qname;
	// Build overlap hash
	std::unordered_map<std::string, BMF::UniqueObservation> obs;
	std::unordered_map<std::string, BMF::UniqueObservation>::iterator found;
	const int sk = 1;
	std::for_each(plp, plp + n_plp, [&found,&obs,&qname](const bam_pileup1_t& plp) {
		if(plp.is_del || plp.is_refskip) return;
		qname = std::string(bam_get_qname(plp.b));
		if((found = obs.find(qname)) == obs.end())
			obs[qname] = BMF::UniqueObservation(plp);
		else found->second.add_obs(plp);
		return;
	});
	// Build vcfline struct
	BMF::VCFPos vcfline = BMF::VCFPos(obs, tid, pos, obs.size());
	vcfline.to_bcf(ret);
	bcf_write(outfp, hdr, ret);
}

int stack_core(aux_t *aux)
{
	bam_plp_t pileup = bam_plp_init((bam_plp_auto_f)read_bam, (void *)aux);
	const bam_pileup1_t *stack;
	bam_plp_set_maxcnt(pileup, aux->max_depth);
	std::vector<khiter_t> sorted_keys = make_sorted_keys(aux->bed);
	int tid, pos, n_plp;
	bcf1_t *v = bcf_init1();
	for(auto& key: sorted_keys) {
		for(uint64_t i = 0; i < kh_val(aux->bed, key).n; ++i) {
			if(aux->iter) hts_itr_destroy(aux->iter);
			const int start = get_start(kh_val(aux->bed, key).intervals[i]);
			const int stop = get_stop(kh_val(aux->bed, key).intervals[i]);
			const int bamtid = (int)kh_key(aux->bed, key);
			aux->iter = bam_itr_queryi(aux->bam_index, kh_key(aux->bed, key), start, stop);
			while((stack = bam_plp_auto(pileup, &tid, &pos, &n_plp)) != 0) {
				if(pos < start && tid == bamtid) continue;
				if(pos >= stop) break;
				process_pileup(v, stack, n_plp, pos, tid, aux->ofp, aux->vh);
			}
		}
	}
	bcf_destroy(v);
	return 0;
}

int stack_main(int argc, char *argv[]) {
	int c;
	unsigned padding = (unsigned)-1;
	aux_t aux = aux_t();
	char *outvcf = NULL;
	std::string bedpath = "";
	int output_bcf = 0;
	if(argc < 2) stack_usage(EXIT_FAILURE);
	while ((c = getopt(argc, argv, "D:q:r:2:S:d:a:s:m:p:f:b:v:o:O:c:BP?hV")) >= 0) {
		switch (c) {
		case 'B': output_bcf = 1; break;
		case 'a': aux.minFA = atoi(optarg); break;
		case 'c': aux.minCount = atoi(optarg); break;
		case 'D': aux.minDuplex = atoi(optarg); break;
		case 's': aux.minFM = atoi(optarg); break;
		case 'm': aux.minMQ = atoi(optarg); break;
		case 'v': aux.minPV = atoi(optarg); break;
		case '2': aux.skip_flag &= (~BAM_FSECONDARY); break;
		case 'S': aux.skip_flag &= (~BAM_FSUPPLEMENTARY); break;
		case 'q': aux.skip_flag &= (~BAM_FQCFAIL); break;
		case 'r': aux.skip_flag &= (~BAM_FDUP); break;
		case 'P': aux.skip_improper = 1; break;
		case 'p': padding = atoi(optarg); break;
		case 'd': aux.max_depth = atoi(optarg); break;
		case 'f': aux.minFR = (float)atof(optarg); break;
		case 'b': bedpath = optarg; break;
		case 'o': outvcf = optarg; break;
		case 'O': aux.minOverlap = atoi(optarg); break;
		//case 'V': aux.vet_all = 1; break;
		case 'h': case '?': stack_usage(EXIT_SUCCESS);
		}
	}
	if(optind >= argc) LOG_ERROR("Insufficient arguments. Input bam required!\n");
	if(padding < 0) {
		LOG_WARNING("Padding not set. Using default %i.\n", DEFAULT_PADDING);
		padding = DEFAULT_PADDING;
	}
	aux.fp = sam_open(argv[optind], "r");
	aux.bh = sam_hdr_read(aux.fp);
	if(!aux.fp || !aux.bh) {
		LOG_ERROR("Could not open input files. Abort!\n");
	}
	aux.bam_index = bam_index_load(aux.fp->fn);
	// TODO: Make BCF header
	aux.vh = bcf_hdr_init(output_bcf ? "wb": "w");
	aux.bed = *bedpath.c_str() ? parse_bed_hash(bedpath.c_str(), aux.bh, padding): build_ref_hash(aux.bh);
	for(auto line: vcf_header_lines)
		if(bcf_hdr_append(aux.vh, line))
			LOG_ERROR("Could not add line %s to header. Abort!\n", line);
	// Add lines to the header for the bed file?
	if(!outvcf) {
		outvcf = (char *)"-";
	}
	aux.ofp = vcf_open(outvcf, output_bcf ? "wb": "w");
	bcf_hdr_write(aux.ofp, aux.vh);
	return stack_core(&aux);
}
