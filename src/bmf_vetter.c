#include "bmf_vetter.h"

int max_depth = 20000;
int padding = 50;


void vetter_error(char *message, int retcode)
{
	fprintf(stderr, message);
	exit(retcode);
}

int vet_func(void *data, bam1_t *b) {
	return 1;
}

/*
 * #define VETTER_OPTIONS \
    {"min-family-agreed",         required_argument, NULL, 'a'}, \
    {"min-family-size",         required_argument, NULL, 's'}, \
    {"min-fraction-agreed",         required_argument, NULL, 'f'}, \
    {"min-mapping-quality",         required_argument, NULL, 'm'}, \
    {"min-phred-quality",         required_argument, NULL, 'p'}, \
    {"in-vcf",         required_argument, NULL, 'v'}, \
    {"out-vcf",         required_argument, NULL, 'o'}, \
    {"bedpath",         required_argument, NULL, 'b'}, \
    {"ref",         required_argument, NULL, 'r'}, \
	{0, 0, 0, 0}
 */

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

void vs_open(vetter_settings_t *settings) {
	htsFormat open_fmt;
	memset(&open_fmt, 0, sizeof(htsFormat));
	open_fmt.category = sequence_data;
	open_fmt.format = bam;
	open_fmt.version.major = 1;
	open_fmt.version.minor = 3;
	vetplp_conf_t *conf = &settings->conf;
	conf->vin = vcf_open(settings->in_vcf_path, settings->vcf_rmode);
	conf->vout = vcf_open(settings->out_vcf_path, settings->vcf_wmode);
	if(strcmp(strrchr(settings->in_vcf_path, '.') - 4, ".vcf.gz") == 0) { // - 4 to skip left for .vcf
		conf->tbx = tbx_index_load(settings->in_vcf_path);
		if(!conf->tbx) {
			fprintf(stderr, "[E:%s] Failed to load index fail for variant file '%s'.", __func__, settings->in_vcf_path);
			exit(EXIT_FAILURE);
		}
#if !NDEBUG
		fprintf(stderr, "[D:%s] Loaded tabix index for bgzipped vcf '%s'.", __func__, settings->in_vcf_path);
#endif
	}
	else if(strcmp(strrchr(settings->in_vcf_path, '.'), ".bcf") == 0) {
		conf->bi = bcf_index_load(settings->in_vcf_path);
#if !NDEBUG
		fprintf(stderr, "[D:%s] Loaded index for bcf file '%s'.", __func__, settings->in_vcf_path);
#endif
	}
	else
		fprintf(stderr, "[W:%s] Not a bcf file or tabixed vcf."
				"Will iterate through the whole file since there's no available index.\n", __func__);

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
	settings->conf.n_iterators = get_nregions(conf->bed);
	if(conf->bed) {
		hts_idx_t *idx = sam_index_load(conf->bam, conf->bam->fn);
		if(!idx) {
			fprintf(stderr, "[%s] Failed to load bam index for file %s. Abort!\n", __func__, conf->bam->fn);
			exit(EXIT_FAILURE);
		}
		int i = 0;
		for(khiter_t k = kh_begin(conf->bed); k != kh_end(conf->bed); ++k) {
			if(!kh_exist(conf->bed, k))
				continue;
			for(int j = 0; j < kh_val(conf->bed, k).n; ++j) {
				// Contig, start, stop
				conf->iterators[i++] = sam_itr_queryi(idx, kh_key(conf->bed, k),
													get_start(kh_val(conf->bed, k).intervals[j]),
													get_stop(kh_val(conf->bed, k).intervals[j]));
			}
		}
		if(i != settings->conf.n_iterators) {
			fprintf(stderr, "The number of bed regions didn't add up. Counted: %i. Expected: %i.", i, settings->conf.n_iterators);
			exit(EXIT_FAILURE);
		}
	}
}

void conf_destroy(vetplp_conf_t *conf)
{
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
	return;
}

void vs_destroy(vetter_settings_t *settings) {
	conf_destroy(&settings->conf);
	free(settings), settings = NULL;
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
		bam_plp_t p = bam_plp_init(&vet_func, (void *)conf);
		hts_itr_t *iter = sam_itr_queryi(conf->bi, rec->rid, rec->pos, rec->pos + 1);
		int ret;
		bam1_t *b = bam_init1();
		while((ret = sam_itr_next(conf->bam, iter, b) >= 0)) {

		}
		hts_itr_destroy(iter);
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
		region_set_t set = kh_val(conf->bed, k);
		const int32_t key = kh_key(conf->bed, k);
		for(uint64_t i = 0; i < set.n; ++i) {
			const uint64_t ivl = set.intervals[i];
			hts_itr_t *bcf_iter = bcf_itr_queryi(conf->vi, key, get_start(ivl), get_stop(ivl) + 1);
			while(bcf_itr_next(conf->vin, bcf_iter, rec) >= 0) {
				if(!bcf_is_snp(rec)) {
					vcf_write(conf->vout, (const bcf_hdr_t *)conf->bh, rec);
					continue;
				}
				bam_plp_t p = bam_plp_init(&vet_func, (void *)conf);
				int ret;
				bam1_t *b = bam_init1();
				bam_destroy1(b);
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
			const uint64_t ivl = cset.intervals[i];
			hts_itr_t *iter = bcf_itr_queryi(conf->vi, key, get_start(ivl), get_start(ivl) + 1);
			// This gets all records in the specified bed region
			while(bcf_itr_next(conf->vin, iter, rec) >= 0) {
				if((conf->bed && !vcf_bed_test(rec, conf->bed)) || !bcf_is_snp(rec))
					continue;
				hts_itr_t *bam_iter = sam_itr_queryi(conf->bi, key, rec->pos, rec->pos + 1);
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
	bam_plp_t iter = bam_plp_init(settings->conf.func, (void *)&settings->conf.params);
	bam_plp_init_overlaps(iter); // Create overlap hashmap for overlapping pairs
	bam_plp_set_maxcnt(iter, max_depth);

	conf->bi = sam_index_load(conf->bam, conf->bam->fn);
	if(!conf->bi) {
		fprintf(stderr, "[%s] Failed to load bam index for file %s. Abort!\n", __func__, conf->bam->fn);
		exit(EXIT_FAILURE);
	}

	if(strcmp(strrchr(settings->conf.vin->fn, '.'), ".vcf") == 0)
		full_iter_loop(settings);
	else
		hts_loop(settings);
	// Main loop
	// Clean up
	bam_plp_destroy(iter);
	return 0;
}

int bmf_vetter_bookends(char *invcf, char *inbam, char *outvcf, char *bed,
						const char *bam_rmode, const char *vcf_rmode,
						const char *vcf_wmode, vparams_t *params)
{
	// Initialize settings struct
	vetter_settings_t *settings = (vetter_settings_t *)calloc(1, sizeof(vetter_settings_t));
	// Initialize
	memcpy(&settings->conf.params, params, sizeof(vparams_t));

	// Copy filenames over and open vcfs.
	strcpy(settings->in_vcf_path, invcf);
	strcpy(settings->bam_path, inbam);
	strcpy(settings->out_vcf_path, outvcf);
	strcpy(settings->bed_path, bed);
	strcpy(settings->vcf_wmode, vcf_wmode);
	strcpy(settings->vcf_rmode, vcf_rmode);
	strcpy(settings->bam_rmode, bam_rmode);

	// Open handles
	vs_open(settings);

	vs_reg_core(settings);
	// Clean up
	vs_destroy(settings);
	free(settings);
	return 0;
}

int bmf_vetter_main(int argc, char *argv[])
{
	const struct option lopts[] = {VETTER_OPTIONS};
	char bam_rmode[3] = "rb";
	char vcf_rmode[4] = "";
	char vcf_wmode[4] = "w";
	char invcf[200] = "";
	char outvcf[200] = "-";
	char inbam[200] = "";
	char bed[200] = "";
	vparams_t params = {
			.minFA = 0u,
			.minPV = 0u,
			.minFM = 0u,
			.minMQ = 0u,
			.minFR = 0.
	};
	int c;
	while ((c = getopt_long(argc, argv, "d:a:s:m:p:f:b:v:o:?h", lopts, NULL)) >= 0) {
		switch (c) {
		case 'a': params.minFA = strtoul(optarg, NULL, 0); break;
		case 's': params.minFM = strtoul(optarg, NULL, 0); break;
		case 'm': params.minMQ = strtoul(optarg, NULL, 0); break;
		case 'v': params.minPV = strtoul(optarg, NULL, 0); break;
		case 'p': padding = atoi(optarg); break;
		case 'd': max_depth = atoi(optarg); break;
		case 'f': params.minFR = atof(optarg); break;
		case 'b': strcpy(bed, optarg); break;
		case 'o': strcpy(outvcf, optarg); break;
		case 'h': /* fall-through */
		case '?': vetter_usage(EXIT_SUCCESS);
		default: vetter_error("Unrecognized option. Abort!\n", EXIT_FAILURE);
		}
	}

	if(argc < 3) {
		fprintf(stderr, "Insufficient arguments. Abort!\n");
		vetter_usage(EXIT_FAILURE);
	}
	if(outvcf[0] == '-') {
		fprintf(stderr, "[%s] Emitting to stdout as vcf.\n", __func__);
		strcpy(vcf_wmode, "w");
	}
	else if(strrchr(outvcf, '.') && strcmp(strrchr(outvcf, '.'), ".bcf") == 0 &&
			!*vcf_rmode)
		strcpy(vcf_rmode, "rb");
	if(!vcf_rmode[0])
		strcpy(vcf_rmode, "r");
	if(!bed[0])
		vetter_error("Bed file required.\n", EXIT_FAILURE);
	strcpy(vcf_wmode, strrchr(outvcf, '.') && strcmp(strrchr(outvcf, '.'), ".bcf") ?
			"w": "wb");
	if(optind + 1 >= argc) {
		vetter_error("Insufficient arguments. Input bam required!\n", EXIT_FAILURE);
	}
	strcpy(invcf, argv[optind]);
	strcpy(inbam, argv[optind + 1]);
	bmf_vetter_bookends(invcf, inbam, outvcf, bed, bam_rmode, vcf_rmode, vcf_wmode, &params);
	return 0;
}
