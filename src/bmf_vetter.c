#include "bmf_vetter.h"

int max_depth = 20000;
int padding = 50;


void vetter_error(char *message, int retcode)
{
	fprintf(stderr, message);
	exit(retcode);
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
	conf->vin = vcf_open(settings->in_vcf_path, "r");
	conf->vout = vcf_open(settings->out_vcf_path, settings->vcf_wmode);
	if(is_bgzipped_vcf(settings->in_vcf_path)) {
		//
		conf->tbx = tbx_index_load(settings->in_vcf_path);
		if(!conf->tbx) {
			LOG_ERROR("Failed to load index fail for variant file '%s'.", settings->in_vcf_path);
		}
		LOG_DEBUG("Loaded tabix index for variant file '%s'.", settings->in_vcf_path);
	}
	else if(file_has_ext(settings->in_vcf_path, "bcf")) {
		conf->bi = bcf_index_load(settings->in_vcf_path);
		LOG_DEBUG("Loaded index for bcf file '%s'.", settings->in_vcf_path);
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

	fprintf(stderr, "Reading header from %s and putting it into a header.\n", settings->conf.bam->fn);
	if((settings->conf.bh = sam_hdr_read(settings->conf.bam)) == NULL)
		vetter_error("Could not read header from bam. Abort!\n", EXIT_FAILURE);
	fprintf(stderr, "Num targets: %.i\n", settings->conf.bh->n_targets);
	conf->vh = bcf_hdr_read(conf->vin);
	bcf_hdr_write(conf->vout, conf->vh);
	if((conf->fai = fai_load(settings->ref_path)) == NULL)
		vetter_error("Could not read Fasta index Abort!\n", EXIT_FAILURE);
	conf->bed = parse_bed_hash(settings->bed_path, conf->bh, padding);
	settings->conf.func = &vet_func;
	settings->conf.n_regions = get_nregions(conf->bed);
}

void conf_destroy(vetplp_conf_t *conf)
{
	if(conf->contig) free(conf->contig), conf->contig = NULL;
	bam_hdr_destroy(conf->bh);
	if(conf->bi)
		hts_idx_destroy(conf->bi);
	if(hts_close(conf->bam))
		vetter_error("Could not close input bam. ??? Abort!\n", EXIT_FAILURE);
	bed_destroy_hash(conf->bed);
	if(hts_close(conf->vin))
		vetter_error("Could not close input vcf. ??? Abort!\n", EXIT_FAILURE);
	if(hts_close(conf->vout))
		vetter_error("Could not close output vcf. ??? Abort!\n", EXIT_FAILURE);
	bcf_hdr_destroy(conf->vh);
	fai_destroy(conf->fai);
	if(conf->tbx) tbx_destroy(conf->tbx), conf->tbx = NULL;
	if(conf->vi) hts_idx_destroy(conf->vi), conf->vi = NULL;
	if(conf->bam_iter) hts_itr_destroy(conf->bam_iter), conf->bam_iter = NULL;
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
		conf->bam_iter = sam_itr_queryi(conf->bi, rec->rid, rec->pos, rec->pos + 1);
		int ret;
		bam1_t *b = bam_init1();
		while((ret = sam_itr_next(conf->bam, conf->bam_iter, b) >= 0)) {

		}
		hts_itr_destroy(conf->bam_iter);
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
	if((conf->vi = bcf_index_load(conf->vin->fn)) == NULL) {
		fprintf(stderr, "Failed to open bcf index. WTF?\n");
		exit(EXIT_FAILURE);
	}
	bcf1_t *rec = bcf_init1();
	bam1_t *b = bam_init1();
	// Each of these iterations sets up a scan for a contig
	for(khint_t ki = kh_begin(conf->bed); ki != kh_end(conf->bed); ++ki) {
		if(ki == kh_end(conf->bed))
			continue;
		const region_set_t cset = kh_val(conf->bed, ki); // Contig set
		const int32_t key = kh_key(conf->bed, ki);
		// Each of these iterations goes through an interval on a contig
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

	}
	bam_destroy1(b);
	bcf_destroy1(rec);
	hts_idx_destroy(conf->vi);
}

int vs_reg_core(vetter_settings_t *settings)
{
	// Set-up
	// Pileup iterator
	vetplp_conf_t *conf = &settings->conf;
	bam_plp_t iter = bam_plp_maxcnt_init(settings->conf.func, (void *)&settings, max_depth);
	bam_plp_init_overlaps(iter); // Create overlap hashmap for overlapping pairs
	conf->pileup = &iter;

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
	return EXIT_SUCCESS;
}

int bmf_vetter_bookends(char *invcf, char *inbam, char *outvcf, char *bed,
						const char *vcf_wmode, vparams_t *params)
{
	int ret;
	// Initialize settings struct
	vetter_settings_t *settings = (vetter_settings_t *)calloc(1, sizeof(vetter_settings_t));
	// Initialize

	settings->conf.minFA = params->minFA;
	settings->conf.minFM = params->minFM;
	settings->conf.minPV = params->minPV;
	settings->conf.minMQ = params->minMQ;
	settings->conf.minFR = params->minFR;
	settings->conf.last_ref_tid = -1; // Make sure it knows it hasn't loaded any yet.
	settings->conf.flag = params->flag;

	// Copy filenames over and open vcfs.
	strcpy(settings->in_vcf_path, invcf);
	strcpy(settings->bam_path, inbam);
	strcpy(settings->out_vcf_path, outvcf);
	strcpy(settings->bed_path, bed);
	strcpy(settings->vcf_wmode, vcf_wmode);

	// Open handles
	vs_open(settings);

	// Run analysis
	ret = vs_reg_core(settings);

	// Clean up
	vs_destroy(settings);
	return ret;
}

int bmf_vetter_main(int argc, char *argv[])
{
	const struct option lopts[] = {VETTER_OPTIONS};
	char vcf_wmode[4] = "w";
	char *invcf = NULL, *outvcf = NULL, *bed = NULL, *inbam = NULL;
	vparams_t params = {
			.minFA = 0u,
			.minPV = 0u,
			.minFM = 0u,
			.minMQ = 0u,
			.minFR = 0.,
			.flag = (BAM_FSECONDARY | BAM_FSUPPLEMENTARY | BAM_FQCFAIL | BAM_FDUP)
	};
	int c;
	while ((c = getopt_long(argc, argv, "q:r:2:$:d:a:s:m:p:f:b:v:o:?h", lopts, NULL)) >= 0) {
		switch (c) {
		case 'a': params.minFA = strtoul(optarg, NULL, 0); break;
		case 's': params.minFM = strtoul(optarg, NULL, 0); break;
		case 'm': params.minMQ = strtoul(optarg, NULL, 0); break;
		case 'v': params.minPV = strtoul(optarg, NULL, 0); break;
		case '2': params.flag &= (~BAM_FSECONDARY); break;
		case '$': params.flag &= (~BAM_FSUPPLEMENTARY); break;
		case 'q': params.flag &= (~BAM_FQCFAIL); break;
		case 'r': params.flag &= (~BAM_FDUP); break;
		case 'p': padding = atoi(optarg); break;
		case 'd': max_depth = atoi(optarg); break;
		case 'f': params.minFR = atof(optarg); break;
		case 'b': strcpy(bed, optarg); break;
		case 'o': strcpy(outvcf, optarg); break;
		case 'h': case '?': vetter_usage(EXIT_SUCCESS);
		}
	}

	if(optind + 1 >= argc)
		vetter_error("Insufficient arguments. Input bam required!\n", EXIT_FAILURE);
	if(outvcf[0] == '-') {
		fprintf(stderr, "[%s] Emitting to stdout as vcf.\n", __func__);
		strcpy(vcf_wmode, "w");
	}
	if(!bed || !*bed)
		vetter_error("Bed file required.\n", EXIT_FAILURE);
	if(file_has_ext(outvcf, "bcf"))
		strcpy(vcf_wmode, "wb");
	invcf = strdup(argv[optind]);
	inbam = strdup(argv[optind + 1]);
	int ret = bmf_vetter_bookends(invcf, inbam, outvcf, bed, vcf_wmode, &params);
	if(ret) {
		fprintf(stderr, "[E:%s:%d] bmf_vetter_bookends returned non-zero exit status '%i'. Abort!\n",
				__func__, __LINE__, ret);
	}
	cond_free(invcf);
	cond_free(outvcf);
	cond_free(bed);
	cond_free(inbam);
	return ret;
}
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
