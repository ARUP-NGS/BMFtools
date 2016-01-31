#include "bmf_vetter.h"

int max_depth = 20000;


void vetter_error(char *message, int retcode)
{
	fprintf(stderr, message);
	exit(retcode);
}


static int read_bam(void *data, bam1_t *b)
{
	aux_t *aux = (aux_t*)data; // data in fact is a pointer to an auxiliary structure
	int ret;
	for(;;)
	{
		ret = aux->iter? sam_itr_next(aux->fp, aux->iter, b) : sam_read1(aux->fp, aux->header, b);
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

/*
 * #define VETTER_OPTIONS \
	{"min-family-agreed",		 required_argument, NULL, 'a'}, \
	{"min-family-size",		 required_argument, NULL, 's'}, \
	{"min-fraction-agreed",		 required_argument, NULL, 'f'}, \
	{"min-mapping-quality",		 required_argument, NULL, 'm'}, \
	{"min-phred-quality",		 required_argument, NULL, 'p'}, \
	{"in-vcf",		 required_argument, NULL, 'v'}, \
	{"out-vcf",		 required_argument, NULL, 'o'}, \
	{"bedpath",		 required_argument, NULL, 'b'}, \
	{"ref",		 required_argument, NULL, 'r'}, \
	{0, 0, 0, 0}
 */
int file_has_ext(char *fn, const char *ext);
int is_bgzipped_vcf(char *fn);
void vetter_usage(int retcode)
{
	char buf[2000];
	sprintf(buf, "Usage:\nbmftools vet <-r/--ref> <ref_path> -o <out.vcf [stdout]> <in.vcf> <in.srt.indexed.bam>\n"
			 "Optional arguments:\n"
			 "-b, --bedpath\tPath to bed file to only validate variants in said region\n"
			 "-s, --min-family-size\tMinimum number of reads in a family to include a that collapsed observation\n"
			 "-f, --min-fraction-agreed\tMinimum fraction of reads in a family agreed on a base call\n"
			 "-v, --min-phred-quality\tMinimum calculated p-value on a base call in phred space\n"
			 "-p, --padding\tNumber of bases outside of bed region to pad.\n"
			 "-a, --min-family-agreed\tMinimum number of reads in a family agreed on a base call\n"
			 "-m, --min-mapping-quality\tMinimum mapping quality for reads for inclusion\n"
			 "Note: fasta reference must be faidx'd.\n");
	vetter_error(buf, retcode);
}

enum vcf_access_mode{
	VCF_UN,
	VCF_BGZIP,
	BCF
};

void vs_open(vetter_settings_t *settings) {
	htsFormat open_fmt = (htsFormat){sequence_data, bam, {1, 3}, gzip, 0, NULL};
	vetplp_conf_t *conf = &settings->conf;
	conf->vin = vcf_open(settings->invcf, "r");
	conf->vout = vcf_open(settings->outvcf, settings->vcf_wmode);
	/*
	if(is_bgzipped_vcf(settings->invcf)) {
		//
		conf->tbx = tbx_index_load(settings->invcf);
		if(!conf->tbx) {
			LOG_ERROR("Failed to load index fail for variant file '%s'.", settings->invcf);
		}
		LOG_DEBUG("Loaded tabix index for variant file '%s'.", settings->invcf);
	}
	else if(file_has_ext(settings->invcf, "bcf")) {
		conf->bi = bcf_index_load(settings->invcf);
		LOG_DEBUG("Loaded index for bcf file '%s'.", settings->invcf);
	} else {
		LOG_WARNING("Not a bcf file or tabixed vcf. "
				"Will iterate through the whole file since there's no available index.\n");
	}

	// Handle bam reading format
	conf->bam = sam_open_format(settings->bam_path, "r", &open_fmt);
	if(conf->bam == NULL) {
		fprintf(stderr, "Failed to open input file. Abort mission!");
		exit(EXIT_FAILURE);
	}
	*/

	fprintf(stderr, "Reading header from %s and putting it into a header.\n", settings->conf.bam->fn);
	if((settings->conf.bh = sam_hdr_read(settings->conf.bam)) == NULL)
		vetter_error("Could not read header from bam. Abort!\n", EXIT_FAILURE);
	fprintf(stderr, "Num targets: %.i\n", settings->conf.bh->n_targets);
	conf->vh = bcf_hdr_read(conf->vin);
	bcf_hdr_write(conf->vout, conf->vh);
}

void conf_destroy(vetplp_conf_t *conf)
{
	//if(conf->contig) free(conf->contig), conf->contig = NULL;
	bam_hdr_destroy(conf->bh);
	if(hts_close(conf->bam))
		vetter_error("Could not close input bam. ??? Abort!\n", EXIT_FAILURE);
	bed_destroy_hash(conf->bed);
	if(hts_close(conf->vin))
		vetter_error("Could not close input vcf. ??? Abort!\n", EXIT_FAILURE);
	if(hts_close(conf->vout))
		vetter_error("Could not close output vcf. ??? Abort!\n", EXIT_FAILURE);
	bcf_hdr_destroy(conf->vh);
	return;
}

void vs_destroy(vetter_settings_t *settings) {
	conf_destroy(&settings->conf);
	free(settings);
}

/* @function
 * :abstract: Iterates over a full VCF, skips variants outside of a region, and
 * fast-forwards the bam to the location for each variant, and tests whether or not that
 * was a correct call based on the supplemental information available in the
 * BMF tags.
 * :param: settings [arg/vetter_settings_t *]
 */
void full_iter_loop(vetter_settings_t *settings)
{
	vetplp_conf_t *conf = &settings->conf;
	bcf1_t *rec = bcf_init1();

	while(bcf_read(conf->vin, conf->vh, rec) >= 0) {
		// Skip over variants outside of our region
		if(conf->bed && !vcf_bed_test(rec, conf->bed))
			continue;
		//conf->bam_iter = sam_itr_queryi(conf->bi, rec->rid, rec->pos, rec->pos + 1);
		int ret;
		bam1_t *b = bam_init1();
		//while((ret = sam_itr_next(conf->bam, conf->bam_iter, b) >= 0)) {

		//}
		//hts_itr_destroy(conf->bam_iter);
		bam_destroy1(b);
	}
	bcf_destroy1(rec);
}


/* @function
 * :abstract: Iterates over a full VCF, skips variants outside of a region, and
 * fast-forwards the bam to the location for each variant, and tests whether or not that
 * was a correct call based on the supplemental information available in the
 * BMF tags.
 * :param: settings [arg/vetter_settings_t *]
 */
void tbx_loop(vetter_settings_t *settings)
{
	vetplp_conf_t *conf = &settings->conf;
	bcf1_t *rec = bcf_init1();
	for(khiter_t k = 0; k != kh_end(conf->bed); ++k) {
		const region_set_t set = kh_val(conf->bed, k);
		const int32_t key = kh_key(conf->bed, k);
		for(uint64_t i = 0; i < set.n; ++i) {
			const uint64_t ivl = set.intervals[i];
			const int start = get_start(ivl);
			const int stop = get_stop(ivl) + 1;

			/*
			conf->bam_iter = sam_itr_queryi(conf->bi, key,
											(start > BAM_FETCH_BUFFER) ? (start - BAM_FETCH_BUFFER): 0,
											stop + BAM_FETCH_BUFFER);
			hts_itr_t *const bcf_iter = bcf_itr_queryi(conf->vi, key, start, stop);
			while(bcf_itr_next(conf->vin, bcf_iter, rec) >= 0) {
				bcf_unpack(rec, BCF_UN_ALL);
				if(!bcf_is_snp(rec)) {
					vcf_write(conf->vout, (const bcf_hdr_t *)conf->bh, rec);
					continue;
				}
				int n_plp, tid, pos;
				const bam_pileup1_t *stack;
				while((stack = bam_plp_auto(*conf->pileup, &tid, &pos, &n_plp)) != NULL) {
					// Note: tid, pos, n_plp are modified in bam_plp_auto.
					// Actual result
				}
			}
			*/
		}
	}
	bcf_destroy1(rec);
}

/* @function
 * :abstract: Creates an iterator over every tabix region in a bed file,
 * fast-forwards the bam to said location, and tests whether or not that
 * was a correct call based on the supplemental information available in the
 * BMF tags.
 * :param: settings [arg/vetter_settings_t *]
 */
void hts_loop(vetter_settings_t *settings)
{
	vetplp_conf_t *conf = &settings->conf;
	/*
	if((conf->vi = bcf_index_load(conf->vin->fn)) == NULL) {
		fprintf(stderr, "Failed to open bcf index. WTF?\n");
		exit(EXIT_FAILURE);
	}
	*/
	bcf1_t *rec = bcf_init1();
	bam1_t *b = bam_init1();
	// Each of these iterations sets up a scan for a contig
	for(khint_t ki = kh_begin(conf->bed); ki != kh_end(conf->bed); ++ki) {
		if(ki == kh_end(conf->bed))
			continue;
		const region_set_t cset = kh_val(conf->bed, ki); // Contig set
		const int32_t key = kh_key(conf->bed, ki);
		// Each of these iterations goes through an interval on a contig
		/*
		for(uint64_t i = 0; i < cset.n; ++i) {
			hts_itr_t *const iter = bcf_itr_queryi(conf->vi, key, get_start(cset.intervals[i]), get_start(cset.intervals[i]) + 1);
			// This gets all records in the specified bed region
			while(bcf_itr_next(conf->vin, (hts_itr_t *)iter, rec) >= 0) {
				if((conf->bed && !vcf_bed_test(rec, conf->bed)) || !bcf_is_snp(rec))
					continue;
				hts_itr_t *const bam_iter = sam_itr_queryi(conf->bi, key, rec->pos - BAM_FETCH_BUFFER, rec->pos + BAM_FETCH_BUFFER + 1);
				while(bam_itr_next(conf->bam, bam_iter, rec) >= 0) {
					if(rec->pos != b->core.pos) {
						fprintf(stderr, "Positions don't equal?? rec: %i. bam: %i.\n", rec->pos, b->core.pos);
						exit(EXIT_FAILURE);
					}
				}
				hts_itr_destroy(bam_iter);
			}
			hts_itr_destroy(iter);
		}
		*/

	}
	bam_destroy1(b);
	bcf_destroy1(rec);
	//hts_idx_destroy(conf->vi);
}
/*
 * typedef struct {
	samFile *fp;
	hts_itr_t *iter;
	bam_hdr_t *header;
	vcfFile *vcf_fp;
	vcfFile *vcf_ofp;
	bcf_hdr_t *vcf_header;
	float minFR; // Minimum fraction of family members agreed on base
	float minAF; // Minimum aligned fraction
	int minFM;
	int minFA;
	int minPV;
	int minMQ;
	int padding;
	int minDuplex;
	int minOverlap;
	int skip_improper;
	uint32_t skip_flag; // Skip reads with any bits set to true
} aux_t;
 */
/*
 * allele here is an unsigned char from seq_nt16_table, meaning that we have converted the variant
 * allele from a character into a bam_seqi format (4 bits per base)
 */
int bmf_pass_var(bcf1_t *vrec, const bam_pileup1_t *plp, unsigned char allele, aux_t *aux, int n_plp,
				int pos) {
	int duplex = 0, overlap = 0, count = 0, i, khr;
	uint32_t *FA1, *PV1, *FA2, *PV2;
	// Build overlap hash
	khash_t(names) *hash = kh_init(names);
	khiter_t k;
	char *qname;
	bam1_t *b;
	uint8_t *seq, *seq2;
	int s, s2;
	char c;
	for(i = 0; i < n_plp; ++i) {
		if(plp[i].is_del || plp[i].is_refskip || plp[i].b->core.flag & BAM_FQCFAIL) continue;
		b = plp[i].b;
		// Skip any reads failed for FA < minFA or FR < minFR
		qname = bam_get_qname(b);
		if((k = kh_get(names, hash, qname)) == kh_end(hash)) {
			kh_put(names, hash, qname, &khr);
			kh_val(hash, k) = &plp[i];
		} else {
			plp[i].b->core.flag &= ~(BAM_FREAD1 | BAM_FREAD2);
			kh_val(hash, k)->b->core.flag |= (BAM_FREAD1 | BAM_FREAD2);
			PV1 = array_tag(kh_val(hash, k)->b, "PV");
			FA1 = array_tag(kh_val(hash, k)->b, "FA");
			seq = bam_get_seq(kh_val(hash, k)->b);
			s = bam_seqi(seq2, kh_val(hash, k)->qpos);
			PV2 = array_tag(plp[i].b, "PV");
			FA2 = array_tag(plp[i].b, "FA");
			seq2 = bam_get_seq(plp[i].b);
			s2 = bam_seqi(seq, plp[i].qpos);
			if(s == s2) {
				PV1[kh_val(hash, k)->qpos] = agreed_pvalues(PV1[kh_val(hash, k)->qpos], PV2[plp[i].qpos]);
				FA1[kh_val(hash, k)->qpos] = FA1[kh_val(hash, k)->qpos] + FA2[plp[i].qpos];
			} else if(s == HTS_N) {
				set_base(seq, seq_nt16_str[bam_seqi(seq2, plp[i].qpos)], kh_val(hash, k)->qpos);
				PV1[kh_val(hash, k)->qpos] = PV2[plp[i].qpos];
				FA1[kh_val(hash, k)->qpos] = FA2[plp[i].qpos];
			} else if(s2 != HTS_N) {
				// Disagreed, both aren't N: N the base, set agrees and p values to 0!
				n_base(seq, kh_val(hash, k)->qpos); // if s2 == HTS_N, do nothing.
				PV1[kh_val(hash, k)->qpos] = 0u;
				FA1[kh_val(hash, k)->qpos] = 0u;
			}
		}
	}
	for(i = 0; i < n_plp; ++i) {
		if(plp[i].is_del || plp[i].is_refskip || (plp[i].b->core.flag & (BAM_FREAD1 | BAM_FREAD2)) == 0) continue;
		b = plp[i].b;
		FA1 = (uint32_t *)array_tag(b, "FA");
		PV1 = (uint32_t *)array_tag(b, "PV");
		if(FA1[plp[i].qpos] < aux->minFA || (float)FA1[plp[i].qpos] / bam_aux2i(bam_aux_get(b, "FM")) < aux->minFR ||
			PV1[plp[i].qpos] < aux->minPV)
			continue;
		seq = bam_get_seq(b);
		if(bam_seqi(seq, plp[i].qpos) == allele) { // Match!
			++count;
			if(bam_aux2i(bam_aux_get(b, "DR"))) ++duplex;
			if((b->core.flag & (BAM_FREAD1 | BAM_FREAD2)) == (BAM_FREAD1 | BAM_FREAD2))
				++overlap;
		}
	}
	return count >= aux->minCount && duplex >= aux->minDuplex && overlap >= aux->minOverlap;
}

/*
 * TODO: Add new tags
 *     1. Number of duplex reads supporting each variant allele
 *
 */

int vet_core(aux_t *aux) {
	khiter_t ki;
	int i, n_plp;
	const bam_pileup1_t *plp = calloc(1, sizeof(bam_pileup1_t*));
	hts_idx_t *idx = sam_index_load(aux->fp, aux->fp->fn);
	tbx_t *vcf_idx = tbx_index_load(aux->vcf_fp->fn);
	bcf1_t *vrec = bcf_init();
	// Unpack all shared data -- up through INFO, but not including FORMAT
	vrec->max_unpack = BCF_UN_FMT;
	hts_itr_t *vcf_iter;
	int32_t *pass_values = (int32_t *)malloc(sizeof(int32_t) * 5);
	for(ki = kh_begin(aux->bed); ki != kh_end(aux->bed); ++ki) {
		if(!kh_exist(aux->bed, ki)) continue;
		for(unsigned j = 0; j < kh_val(aux->bed, ki).n; ++j) {

			int tid, start, stop, pos;

			// Handle coordinates
			tid = kh_key(aux->bed, ki);
			start = get_start(kh_val(aux->bed, ki).intervals[j]);
			stop = get_stop(kh_val(aux->bed, ki).intervals[j]);

			vcf_iter = tbx_itr_queryi(vcf_idx, tid, start, stop);
			if (aux->iter) hts_itr_destroy(aux->iter);
			aux->iter = sam_itr_queryi(idx, tid, start, stop);
			bam_plp_t pileup = bam_plp_init(read_bam, (void *)aux);
			bam_plp_set_maxcnt(pileup, max_depth);
			while(bcf_itr_next(aux->vcf_fp, vcf_iter, vrec) >= 0) {
				while(vrec->pos < start && bcf_itr_next(aux->vcf_fp, vcf_iter, vrec) >= 0) {
					/* Zoom ahead until you're at the correct position */
				}
				if(!bcf_is_snp(vrec)) continue; // Only handle simple SNVs
				bcf_unpack(vrec, BCF_UN_STR); // Unpack the allele fields
				while (pos < vrec->pos && ((plp = bam_plp_auto(pileup, &tid, &pos, &n_plp)) != 0)) {
					/* Zoom ahead until you're at the correct position */
				}
				if(pos != vrec->pos) {
					LOG_INFO("Position %i (1-based) not found in pileups in bam. Super weird...\n", vrec->pos);
					return -1;
				}
				// Check each variant

				for(unsigned i = 0; i < vrec->n_allele; ++i)
					pass_values[i] = bmf_pass_var(vrec, plp, seq_nt16_table[(uint8_t)*(vrec->d.allele[i])], aux, n_plp, pos);
				bcf_update_format(aux->vcf_header, vrec, "BMF", (void *)pass_values, vrec->n_allele, BCF_HT_INT);

				// Pass or fail them individually.
				bcf_write(aux->vcf_ofp, aux->vcf_header, vrec);
			}
			tbx_itr_destroy(vcf_iter);
			bam_plp_destroy(pileup);
		}
	}
	free(pass_values);
	hts_idx_destroy(idx);
	bcf_destroy(vrec);
	return EXIT_SUCCESS;
}

int vs_reg_core(vetter_settings_t *settings)
{
	/*
	// Set-up
	// Pileup iterator
	vetplp_conf_t *conf = &settings->conf;
	bam_plp_t iter = bam_plp_maxcnt_init(settings->conf.func, (void *)&settings, max_depth);
	bam_plp_init_overlaps(iter); // Create overlap hashmap for overlapping pairs
	//conf->pileup = &iter;

	conf->bi = sam_index_load(conf->bam, conf->bam->fn);
	if(!conf->bi) {
		fprintf(stderr, "[%s] Failed to load bam index for file %s. Abort!\n", __func__, conf->bam->fn);
		exit(EXIT_FAILURE);
	}

	if(strcmp(strrchr(settings->conf.vin->fn, '.'), ".vcf") == 0)
		full_iter_loop(settings);
	else if(strcmp(strrchr(settings->conf.vin->fn, '.'), ".bcf") == 0)
		hts_loop(settings);
	else if(strcmp(strrchr(settings->conf.vin->fn, '.') - 4, ".vcf.gz"))
		tbx_loop(settings);
	// Main loop
	// Clean up
	bam_plp_destroy(iter);
	*/
	return EXIT_SUCCESS;
}

int vetter_main(int argc, char *argv[])
{
	const struct option lopts[] = {VETTER_OPTIONS};
	char vcf_wmode[4] = "w";
	char *outvcf = NULL, *bed = NULL;
	int c;
	int padding = 0;
	htsFormat open_fmt = (htsFormat){sequence_data, bam, {1, 3}, gzip, 0, NULL};
	/* 	float minFR; // Minimum fraction of family members agreed on base
	float minAF; // Minimum aligned fraction
	int minFM;
	int minFA;
	int minPV;
	int minMQ;
	int padding;
	int minDuplex;
	int minOverlap;
	uint32_t skip_flag; // Skip reads with any bits set to true*/
	aux_t aux = {0};
	aux.max_depth = (1 << 18); // Default max depth

	//while ((c = getopt_long(argc, argv, "q:r:2:$:d:a:s:m:p:f:b:v:o:O:c:P?h", lopts, NULL)) >= 0) {
	while ((c = getopt(argc, argv, "q:r:2:$:d:a:s:m:p:f:b:v:o:O:c:P?h")) >= 0) {
		switch (c) {
		case 'a': aux.minFA = atoi(optarg); break;
		case 'c': aux.minCount = atoi(optarg); break;
		case 's': aux.minFM = atoi(optarg); break;
		case 'm': aux.minMQ = atoi(optarg); break;
		case 'v': aux.minPV = atoi(optarg); break;
		case '2': aux.skip_flag &= (~BAM_FSECONDARY); break;
		case '$': aux.skip_flag &= (~BAM_FSUPPLEMENTARY); break;
		case 'q': aux.skip_flag &= (~BAM_FQCFAIL); break;
		case 'r': aux.skip_flag &= (~BAM_FDUP); break;
		case 'P': aux.skip_improper = 1; break;
		case 'p': padding = atoi(optarg); break;
		case 'd': aux.max_depth = atoi(optarg); break;
		case 'f': aux.minFR = (float)atof(optarg); break;
		case 'b': bed = strdup(optarg); break;
		case 'o': outvcf = strdup(optarg); break;
		case 'O': aux.minOverlap = atoi(optarg); break;
		case 'h': case '?': vetter_usage(EXIT_SUCCESS);
		}
	}

	// Check for required tags.
	if(aux.minAF)
		check_bam_tag_exit(argv[optind + 1], "AF");
	check_bam_tag_exit(argv[optind + 1], "FA");
	check_bam_tag_exit(argv[optind + 1], "FM");
	check_bam_tag_exit(argv[optind + 1], "FP");
	check_bam_tag_exit(argv[optind + 1], "PV");
	check_bam_tag_exit(argv[optind + 1], "RV");



	if(optind + 1 >= argc)
		vetter_error("Insufficient arguments. Input bam required!\n", EXIT_FAILURE);
	if(outvcf[0] == '-') {
		LOG_INFO("Emitting to stdout as vcf.\n");
		strcpy(vcf_wmode, "w");
	}
	if(!bed || !*bed)
		vetter_error("Bed file required.\n", EXIT_FAILURE);
	if(file_has_ext(outvcf, "bcf"))
		strcpy(vcf_wmode, "wb");

	// Open bam
	aux.fp = sam_open_format(argv[optind + 1], "r", &open_fmt);
	if(!aux.fp) {LOG_ERROR("Could not open input bam %s. Abort!\n", argv[optind + 1]);}
	aux.header = sam_header_read(aux.fp);

	// Open input vcf
	if(!aux.header || aux.header->n_targets == 0) {LOG_ERROR("Could not read header from bam %s. Abort!\n", argv[optind + 1]);}
	aux.vcf_fp = vcf_open(argv[optind], "r");
	if(!aux.vcf_fp) {LOG_ERROR("Could not open input [bv]cf %s. Abort!\n", argv[optind]);}
	aux.vcf_header = vcf_hdr_read(aux.vcf_fp);
	if(!aux.vcf_header) {LOG_ERROR("Could not read header from input [bv]cf %s. Abort!\n", argv[optind]);}
	// Open bed file
	aux.bed = parse_bed_hash(bed, aux.header, padding);

	// Add lines to header
	size_t n_header_lines = COUNT_OF(bmf_header_lines);
	for(int i = 0; i < n_header_lines; ++i) bcf_hdr_append(aux.vcf_header, bmf_header_lines[i]);
	bcf_hdr_printf(aux.vcf_header, "##bed_filename=\"%s\"", bed);
	{ // New block so tmpstr is cleared
		kstring_t tmpstr = {0};
		ksprintf(&tmpstr, "##cmdline=");
		for(int i = 0; i < argc; ++i) kputs(argv[i], &tmpstr), kputc(' ', &tmpstr);
		bcf_hdr_append(aux.vcf_header, tmpstr.s);
		free(tmpstr.s);
	}
	bcf_hdr_printf(aux.vcf_header, "##bmftools_version=\"%s\"", VERSION);
	bcf_hdr_write(aux.vcf_ofp, aux.vcf_header);

	// Open out vcf
	int ret = vet_core(&aux);
	if(ret)
		fprintf(stderr, "[E:%s:%d] vet_core returned non-zero exit status '%i'. Abort!\n",
				__func__, __LINE__, ret);
	cond_free(outvcf);
	cond_free(bed);
	return ret;
}
/*
static int vet_func(void *data, bam1_t *b)
{
	vetplp_conf_t *conf = (vetplp_conf_t *)data;
	int ret, skip = 0, ref_len;
	do {
		ret = conf->bam_iter ? sam_itr_next(conf->bam, conf->bam_iter, b) : sam_read1(conf->bam, conf->bh, b);
		if (ret < 0) break;
		// The 'B' cigar operation is not part of the specification, considering as obsolete.
		//  bam_remove_B(b);
		if (b->core.tid < 0 || (b->core.flag&(BAM_FUNMAP))) { // exclude unmapped and qc fail reads.
			skip = 1;
			continue;
		}
		if (conf->bed) { // test overlap
			skip = !bed_test(b, conf->bed);
			if (skip) continue;
		}

		if (conf->fai && b->core.tid >= 0) {
			if(conf->last_ref_tid != b->core.tid) {
				if(conf->contig) free(conf->contig);
				conf->contig = fai_fetch(conf->fai, conf->bh->target_name[b->core.tid], &ref_len);
			}
			if (ref_len <= b->core.pos) { // exclude reads outside of the reference sequence
				fprintf(stderr,"[%s] Skipping because %d is outside of %d [ref:%d]\n",
						__func__, b->core.pos, ref_len, b->core.tid);
				skip = 1;
				continue;
			}
		}

		skip = 0;
		if (b->core.qual < conf->minMQ) skip = 1;
		else if((conf->flag & (b->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY | BAM_FQCFAIL | BAM_FDUP))) ||
				((conf->flag & (SKIP_IMPROPER)) && (b->core.flag&BAM_FPAIRED) && b->core.flag&BAM_FPROPER_PAIR) ||
				(bam_aux2i(bam_aux_get(b, "FM")) < conf->minFM)) skip = 1;
	} while (skip);
	return ret;
}
*/
