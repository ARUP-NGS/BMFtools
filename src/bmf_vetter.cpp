#include "bmf_vetter.h"

int max_depth = (1 << 18); // 262144


void vetter_error(const char *message, int retcode)
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
void vetter_usage(int retcode)
{
	char buf[2000];
	sprintf(buf,
			"Usage:\nbmftools vet -o <out.vcf [stdout]> <in.vcf.gz/in.bcf> <in.srt.indexed.bam>\n"
			"Optional arguments:\n"
			"-b, --bedpath\tPath to bed file to only validate variants in said region. REQUIRED.\n"
			"-s, --min-family-size\tMinimum number of reads in a family to include a that collapsed observation\n"
			"-2, --skip-secondary\tSkip secondary alignments.\n"
			"-S, --skip-supplementary\tSkip supplementary alignments.\n"
			"-q, --skip-qcfail\tSkip reads marked as QC fail.\n"
			"-f, --min-fraction-agreed\tMinimum fraction of reads in a family agreed on a base call\n"
			"-v, --min-phred-quality\tMinimum calculated p-value on a base call in phred space\n"
			"-p, --padding\tNumber of bases outside of bed region to pad.\n"
			"-a, --min-family-agreed\tMinimum number of reads in a family agreed on a base call\n"
			"-m, --min-mapping-quality\tMinimum mapping quality for reads for inclusion\n"
			"-B, --emit-bcf-format\tEmit bcf-formatted output. (Defaults to vcf).\n"
			);
	vetter_error(buf, retcode);
}

enum vcf_access_mode{
	VCF_UN,
	VCF_BGZIP,
	BCF
};

/*
 * allele here is an unsigned char from seq_nt16_table, meaning that we have converted the variant
 * allele from a character into a bam_seqi format (4 bits per base)
 */
int bmf_pass_var(bcf1_t *vrec, const bam_pileup1_t *plp, unsigned char allele, aux_t *aux, int n_plp,
				int pos) {
	int duplex = 0, overlap = 0, count = 0, i, khr, s, s2;
	khiter_t k;
	uint32_t *FA1, *PV1, *FA2, *PV2;
	char *qname;
	bam1_t *b;
	uint8_t *seq, *seq2, *tmptag, *drdata;
	// Build overlap hash
	khash_t(names) *hash = kh_init(names);
	const int sk = 1;
	// Set the r1/r2 flags for the reads to ignore to 0
	// Set the ones where we see it twice to (BAM_FREAD1 | BAM_FREAD2).
	for(i = 0; i < n_plp; ++i) {
		if(plp[i].is_del || plp[i].is_refskip) continue;
		b = plp[i].b;
		// Skip any reads failed for FA < minFA or FR < minFR
		qname = bam_get_qname(b);
		k = kh_get(names, hash, qname);
		if(k == kh_end(hash)) {
			kh_put(names, hash, qname, &khr);
			k = kh_get(names, hash, qname);
			kh_val(hash, k) = &plp[i];
		} else {
			bam_aux_append(plp[i].b, "SK", 'i', sizeof(int), (uint8_t *)&sk); // Skip
			bam_aux_append(kh_val(hash, k)->b, "KR", 'i', sizeof(int), (uint8_t *)&sk); // Keep Read
			PV1 = (uint32_t *)array_tag(kh_val(hash, k)->b, "PV");
			FA1 = (uint32_t *)array_tag(kh_val(hash, k)->b, "FA");
			seq = bam_get_seq(kh_val(hash, k)->b);
			s = bam_seqi(seq, kh_val(hash, k)->qpos);
			PV2 = (uint32_t *)array_tag(plp[i].b, "PV");
			FA2 = (uint32_t *)array_tag(plp[i].b, "FA");
			seq2 = bam_get_seq(plp[i].b);
			s2 = bam_seqi(seq2, plp[i].qpos);
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
		if(plp[i].is_del || plp[i].is_refskip) continue;
		if((tmptag = bam_aux_get(plp[i].b, "SK")) != NULL) {
			// If it has the SK tag, get rid of it, but skip in the pileup.
			// That way, each position isn't affecting the results of neighboring calls.
			bam_aux_del(plp[i].b, tmptag);
			continue;
		}
		b = plp[i].b;
		FA1 = (uint32_t *)array_tag(b, "FA");
		PV1 = (uint32_t *)array_tag(b, "PV");
		if(FA1[plp[i].qpos] < aux->minFA || (float)FA1[plp[i].qpos] / bam_aux2i(bam_aux_get(b, "FM")) < aux->minFR ||
			PV1[plp[i].qpos] < aux->minPV)
			continue;
		seq = bam_get_seq(b);
		if(bam_seqi(seq, plp[i].qpos) == allele) { // Match!
			++count;
			if((drdata = bam_aux_get(b, "DR")) != NULL &&
				bam_aux2i(drdata)) ++duplex; // Has DR tag and its value is nonzero.
			if((tmptag = bam_aux_get(b, "KR")) != NULL) {
				++overlap;
				bam_aux_del(b, tmptag);
			}
		}
	}
	kh_destroy(names, hash);
	return count >= aux->minCount && duplex >= aux->minDuplex && overlap >= aux->minOverlap;
}

/*
 * TODO: Add new tags
 *     1. Number of duplex reads supporting each variant allele.
 *     2. Number of passing reads for each allele.
 *     3. Number of unique observations (subtracting overlapping reads) for each allele.
 */

int read_bcf(aux_t *aux, hts_itr_t *vcf_iter, bcf1_t *vrec, int start, int tid)
{
	int ret;
	if(vcf_iter) {
		ret = bcf_itr_next(aux->vcf_fp, vcf_iter, vrec);
		LOG_DEBUG("Accessing vcf_iter at %p. tid: %i. start: %i. vrec->pos: %i.\n", (void *)vcf_iter, tid, start, vrec->pos);
		if((ret = bcf_itr_next(aux->vcf_fp, vcf_iter, vrec)) < 0) return ret;
		LOG_DEBUG("pos: %i.\n", vrec->pos);
		while(vrec->pos < start)
		{
			/* Zoom ahead until you're at the correct position */
			if((ret = bcf_itr_next(aux->vcf_fp, vcf_iter, vrec)) < 0) return ret;
		}
	} else {
		if((ret = bcf_read(aux->vcf_fp, aux->vcf_header, vrec)) < 0) return ret;
		while(vrec->rid < tid) {
			ret = bcf_read(aux->vcf_fp, aux->vcf_header, vrec);
		}
		while(vrec->pos < start) {
			if(vrec->rid != tid) return ret;
			ret = bcf_read(aux->vcf_fp, aux->vcf_header, vrec);
		}
	}
	return ret;
}

int vet_core(aux_t *aux) {
	khiter_t ki;
	int n_plp, vcf_iter_ret;
	const bam_pileup1_t *plp;
	tbx_t *vcf_idx = NULL;
	hts_idx_t *bcf_idx = NULL;
	hts_idx_t *idx = sam_index_load(aux->fp, aux->fp->fn);
	switch(hts_get_format(aux->vcf_fp)->format) {
	case vcf:
		vcf_idx = tbx_index_load(aux->vcf_fp->fn);
		if(!vcf_idx) {
			LOG_WARNING("Could not load TBI index for %s. Defaulting to all variants.\n", aux->vcf_fp->fn);
		}
		break;
	case bcf:
		bcf_idx = bcf_index_load(aux->vcf_fp->fn);
		if(!bcf_idx) {
			LOG_ERROR("Could not load CSI index: %s\n", aux->vcf_fp->fn);
		}
		break;
	default:
		LOG_ERROR("Unrecognized variant file type! (%i).\n", hts_get_format(aux->vcf_fp)->format);
	}
	/*
	if(!(vcf_idx || bcf_idx)) {
		LOG_ERROR("Require an indexed variant file. Abort!\n");
	}
	*/
	bcf1_t *vrec = bcf_init();
	// Unpack all shared data -- up through INFO, but not including FORMAT
	vrec->max_unpack = BCF_UN_FMT;
	hts_itr_t *vcf_iter = NULL;
	int32_t *pass_values = (int32_t *)malloc(sizeof(int32_t) * 5);
	bam_plp_t pileup = bam_plp_init(read_bam, (void *)aux);
	bam_plp_set_maxcnt(pileup, max_depth);
	for(ki = kh_begin(aux->bed); ki != kh_end(aux->bed); ++ki) {
		hash_label:
		if(!kh_exist(aux->bed, ki)) continue;
		for(unsigned j = 0; j < kh_val(aux->bed, ki).n; ++j) {
			int tid, start, stop, pos = -1;

			// Handle coordinates
			tid = kh_key(aux->bed, ki);
			start = get_start(kh_val(aux->bed, ki).intervals[j]);
			stop = get_stop(kh_val(aux->bed, ki).intervals[j]);
			LOG_DEBUG("Beginning to work through region #%i on contig %s:%i-%i.\n", j + 1, aux->header->target_name[tid], start, stop);

			// Fill vcf_iter from tbi or csi index. If both are null, go through the full file.
			vcf_iter = vcf_idx ? hts_itr_query((const hts_idx_t *)vcf_idx, tid, start, stop, &tbx_readrec): bcf_idx ? bcf_itr_queryi(bcf_idx, tid, start, stop): NULL;

			while((vcf_iter_ret = read_bcf(aux, vcf_iter, vrec, start, tid)) >= 0) {
				if(!bcf_is_snp(vrec) || !vcf_bed_test(vrec, aux->bed)) {
					bcf_write(aux->vcf_ofp, aux->vcf_header, vrec);
					continue; // Only handle simple SNVs
				}
				if(vrec->rid > tid) {
					// No variants on this contig.
					// Increment ki to go to the next contig.
					++ki;
					goto hash_label;
				}
				bcf_unpack(vrec, BCF_UN_STR); // Unpack the allele fields
				if (aux->iter) hts_itr_destroy(aux->iter);
				aux->iter = sam_itr_queryi(idx, vrec->rid, vrec->pos, stop);
				plp = bam_plp_auto(pileup, &tid, &pos, &n_plp);
				if(!plp) {
					LOG_ERROR("Could not make pileup for region %s:%i-%i.\n", aux->header->target_name[tid], start, stop);
				}
				while ((tid < vrec->rid || pos < vrec->pos) && ((plp = bam_plp_auto(pileup, &tid, &pos, &n_plp)) != 0)) {
					/* Zoom ahead until you're at the correct position */
				}
				LOG_DEBUG("tid: %i. rid: %i. var pos: %i.\n", tid, vrec->rid, vrec->pos);
				if(pos != vrec->pos || tid != vrec->rid) {
					LOG_WARNING("Position %s:%i (1-based) not found in pileups in bam. Writing unmodified. Super weird...\n", aux->header->target_name[tid], vrec->pos + 1);
					bcf_write(aux->vcf_ofp, aux->vcf_header, vrec);
					continue;
				}
				// Check each variant

				for(unsigned i = 0; i < vrec->n_allele; ++i)
					pass_values[i] = bmf_pass_var(vrec, plp, seq_nt16_table[(uint8_t)(vrec->d.allele[i][0])], aux, n_plp, pos);
				LOG_DEBUG("n allele: %i.\n", vrec->n_allele);
				bcf_update_info(aux->vcf_header, vrec, "BMF_VET", (void *)pass_values, vrec->n_allele, BCF_HT_INT);

				// Pass or fail them individually.
				bcf_write(aux->vcf_ofp, aux->vcf_header, vrec);
			}
			if(vcf_iter) tbx_itr_destroy(vcf_iter);
		}
	}
	bam_plp_destroy(pileup);
	if(bcf_idx) hts_idx_destroy(bcf_idx);
	if(vcf_idx) tbx_destroy(vcf_idx);
	if(aux->iter) hts_itr_destroy(aux->iter);
	free(pass_values);
	hts_idx_destroy(idx);
	bcf_destroy(vrec);
	return EXIT_SUCCESS;
}

int vetter_main(int argc, char *argv[])
{
	const struct option lopts[] = {VETTER_OPTIONS};
	char vcf_wmode[4] = "w";
	char *outvcf = NULL, *bed = NULL;
	int c;
	int padding = 0, output_bcf = 0;
	// Defaults to outputting textual (vcf)
	htsFormat open_fmt = {sequence_data, bam, {1, 3}, gzip, 0, NULL};
	aux_t aux = {0};
	aux.max_depth = (1 << 18); // Default max depth
	if(argc < 3) vetter_usage(EXIT_FAILURE);

	while ((c = getopt_long(argc, argv, "D:q:r:2:S:d:a:s:m:p:f:b:v:o:O:c:BP?h", lopts, NULL)) >= 0) {
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
		case 'b': bed = strdup(optarg); break;
		case 'o': outvcf = strdup(optarg); break;
		case 'O': aux.minOverlap = atoi(optarg); break;
		case 'h': case '?': vetter_usage(EXIT_SUCCESS);
		}
	}

	// Check for required tags.
	if(aux.minAF) check_bam_tag_exit(argv[optind + 1], "AF");
	check_bam_tag_exit(argv[optind + 1], "FA");
	check_bam_tag_exit(argv[optind + 1], "FM");
	check_bam_tag_exit(argv[optind + 1], "FP");
	check_bam_tag_exit(argv[optind + 1], "PV");
	check_bam_tag_exit(argv[optind + 1], "RV");



	if(optind + 1 >= argc)
		vetter_error("Insufficient arguments. Input bam required!\n", EXIT_FAILURE);
	strcpy(vcf_wmode, output_bcf ? "wb": "w");
	if(!outvcf) // Default to emitting to stdout.
		outvcf = strdup("-");
	if(strcmp(outvcf, "-") == 0) {
		LOG_INFO("Emitting to stdout in %s format.\n", output_bcf ? "bcf": "vcf");
	}
	// Open bam
	aux.fp = sam_open_format(argv[optind + 1], "r", &open_fmt);
	if(!aux.fp) {LOG_ERROR("Could not open input bam %s. Abort!\n", argv[optind + 1]);}
	aux.header = sam_header_read(aux.fp);

	// Open input vcf
	if(!aux.header || aux.header->n_targets == 0) {LOG_ERROR("Could not read header from bam %s. Abort!\n", argv[optind + 1]);}
	// Open bed file
	// if no bed provided, do whole genome.
	if(!bed) {
		LOG_ERROR("No bed file provided. Required. Abort!\n");
		//bed = strdup("FullGenomeAnalysis");
		//LOG_WARNING("No bed file provided. Defaulting to whole genome analysis.\n");
		//aux.bed = build_ref_hash(aux.header);
	} else aux.bed = parse_bed_hash(bed, aux.header, padding);
	//check_vcf_open(argv[optind], aux.vcf_fp, aux.vcf_header);
	aux.vcf_fp = vcf_open(argv[optind], "r");
	if(!aux.vcf_fp) {
		LOG_ERROR("BLAH");
	}
	aux.vcf_header = bcf_hdr_read(aux.vcf_fp);
	if(!aux.vcf_header) {
		LOG_ERROR("BLAH");
	}

	// Add lines to header
	for(unsigned i = 0; i < COUNT_OF(bmf_header_lines); ++i)
		if(bcf_hdr_append(aux.vcf_header, bmf_header_lines[i]))
			fprintf(stderr, "[E:%s:%d] Could not add header line '%s'. Abort!\n", __func__, __LINE__, bmf_header_lines[i]), exit(EXIT_FAILURE);
	bcf_hdr_printf(aux.vcf_header, "##bed_filename=\"%s\"", bed);
	{ // New block so tmpstr is cleared
		kstring_t tmpstr = {0};
		ksprintf(&tmpstr, "##cmdline=");
		kputs("bmftools ", &tmpstr);
		for(int i = 0; i < argc; ++i) kputs(argv[i], &tmpstr), kputc(' ', &tmpstr);
		bcf_hdr_append(aux.vcf_header, tmpstr.s);
		free(tmpstr.s);
	}
	bcf_hdr_printf(aux.vcf_header, "##bmftools_version=\"%s\"", VERSION);

	// Open output vcf
	aux.vcf_ofp = vcf_open(outvcf, vcf_wmode);
	if(!aux.vcf_ofp) {
		LOG_ERROR("Could not open output vcf '%s' for writing. Abort!\n", outvcf);
	}
	bcf_hdr_write(aux.vcf_ofp, aux.vcf_header);

	// Open out vcf
	int ret = vet_core(&aux);
	if(ret)
		fprintf(stderr, "[E:%s:%d] vet_core returned non-zero exit status '%i'. Abort!\n",
				__func__, __LINE__, ret);
	sam_close(aux.fp);
	bam_hdr_destroy(aux.header);
	vcf_close(aux.vcf_fp);
	vcf_close(aux.vcf_ofp);
	bcf_hdr_destroy(aux.vcf_header);
	bed_destroy_hash(aux.bed);
	cond_free(outvcf);
	cond_free(bed);
	LOG_INFO("Successfully completed!\n");
	return ret;
}
