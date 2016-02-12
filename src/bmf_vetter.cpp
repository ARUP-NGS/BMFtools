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

void vetter_usage(int retcode)
{
	const char *buf =
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
			"-B, --emit-bcf-format\tEmit bcf-formatted output. (Defaults to vcf).\n";
	vetter_error(buf, retcode);
}

inline int arr_qpos(const bam_pileup1_t *plp)
{
	LOG_DEBUG("qpos: %i.\n", plp->qpos);
	LOG_DEBUG("l_qseq: %i.\n", plp->b->core.l_qseq);
	LOG_DEBUG("Arr qpos: %i.\n", (plp->b->core.flag & BAM_FREVERSE) ? plp->b->core.l_qseq - 1 - plp->qpos: plp->qpos);
	return (plp->b->core.flag & BAM_FREVERSE) ? plp->b->core.l_qseq - 1 - plp->qpos: plp->qpos;
}

/*
 * :param: [bcf1_t *] vrec - Variant record to test.
 *   # UniObs passing
 * 2. What do I want for INFO?
 *
 */
void bmf_var_tests(bcf1_t *vrec, const bam_pileup1_t *plp, int n_plp, aux_t *aux, std::vector<int>& pass_values,
		std::vector<int>& n_obs, std::vector<int>& n_duplex, std::vector<int>& n_overlaps, std::vector<int> &n_failed,
		int& n_all_overlaps, int& n_all_duplex, int& n_all_disagreed) {
	int khr, s, s2;
	n_all_disagreed = n_all_overlaps = 0;
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
			++n_all_overlaps;
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
				++n_all_disagreed;
				// Disagreed, both aren't N: N the base, set agrees and p values to 0!
				n_base(seq, kh_val(hash, k)->qpos); // if s2 == HTS_N, do nothing.
				PV1[arr_qpos1] = 0u;
				FA1[arr_qpos1] = 0u;
			}
		}
	}
	for(int j = 0; j < vrec->n_allele; ++j) {
		for(int i = 0; i < n_plp; ++i) {
			if(plp[i].is_del || plp[i].is_refskip) continue;
			if((tmptag = bam_aux_get(plp[i].b, "SK")) != NULL) {
				continue;
			}

			seq = bam_get_seq(plp[i].b);
			FA1 = (uint32_t *)array_tag(plp[i].b, "FA");
			PV1 = (uint32_t *)array_tag(plp[i].b, "PV");
#if !NDEBUG
			fprintf(stderr, "Read name: %s.\n", bam_get_qname(plp[i].b));
			for(int k1 = 0; k1 < plp[i].b->core.l_qseq; ++k1) {
				fprintf(stderr, ",%u", PV1[k1]);
			}
			fputc('\n', stderr);
			for(int k1 = 0; k1 < plp[i].b->core.l_qseq; ++k1) {
				fprintf(stderr, ",%u", FA1[k1]);
			}
			fputc('\n', stderr);
#endif
			if(bam_seqi(seq, plp[i].qpos) == seq_nt16_table[(uint8_t)vrec->d.allele[j][0]]) { // Match!
				const int32_t arr_qpos1 = arr_qpos(&plp[i]);
				if(FA1[arr_qpos1] < aux->minFA || PV1[arr_qpos1] < aux->minPV ||
						(float)FA1[arr_qpos1] / bam_aux2i(bam_aux_get(plp[i].b, "FM")) < aux->minFR) {
					++n_failed[j];
					continue;
				}
				LOG_DEBUG("Note: PV1[%i] value (%u) has to be greater than minPV now. (%u)\n", arr_qpos1, PV1[arr_qpos1], aux->minPV);
				++n_obs[j];
				if((drdata = bam_aux_get(plp[i].b, "DR")) != NULL && bam_aux2i(drdata)) {
					++n_duplex[j]; // Has DR tag and its value is nonzero.
				}
				if((tmptag = bam_aux_get(plp[i].b, "KR")) != NULL) {
					++n_overlaps[j];
					bam_aux_del(plp[i].b, tmptag);
				}
			}
		}
		pass_values[j] = n_obs[j] >= aux->minCount && n_duplex[j] >= aux->minDuplex && n_overlaps[j] >= aux->minOverlap;

	}
	for(int i = 0; i < n_plp; ++i) {
		if((tmptag = bam_aux_get(plp[i].b, "SK")) != NULL) bam_aux_del(plp[i].b, tmptag);
	}
	kh_destroy(names, hash);
	n_all_duplex = std::accumulate(n_duplex.begin(), n_duplex.begin() + vrec->n_allele, 0);
	//return count >= aux->minCount && duplex >= aux->minDuplex && overlap >= aux->minOverlap;
}

/*
 * allele here is an unsigned char from seq_nt16_table, meaning that we have converted the variant
 * allele from a character into a bam_seqi format (4 bits per base)
 */
int bmf_pass_var(bcf1_t *vrec, const bam_pileup1_t *plp, unsigned char allele, aux_t *aux, int n_plp,
				int pos, int *n_disagreed, int *n_overlap, int *n_duplex, int *n_uniobs) {
	int duplex = 0, overlap = 0, count = 0, i, khr, s, s2;
	khiter_t k;
	uint32_t *FA1, *PV1, *FA2, *PV2;
	char *qname;
	uint8_t *seq, *seq2, *tmptag, *drdata;
	*n_overlap = *n_disagreed = 0;
	// Build overlap hash
	khash_t(names) *hash = kh_init(names);
	const int sk = 1;
	// Set the r1/r2 flags for the reads to ignore to 0
	// Set the ones where we see it twice to (BAM_FREAD1 | BAM_FREAD2).
	for(i = 0; i < n_plp; ++i) {
		if(plp[i].is_del || plp[i].is_refskip) continue;
		// Skip any reads failed for FA < minFA or FR < minFR
		qname = bam_get_qname(plp[i].b);
		k = kh_get(names, hash, qname);
		if(k == kh_end(hash)) {
			k = kh_put(names, hash, qname, &khr);
			kh_val(hash, k) = &plp[i];
		} else {
			const int32_t arr_qpos1 = arr_qpos(kh_val(hash, k));
			const int32_t arr_qpos2 = arr_qpos(&plp[i]);
			++*n_overlap;
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
				PV1[arr_qpos1] = agreed_pvalues(PV1[arr_qpos1], PV2[arr_qpos2]);
				FA1[arr_qpos1] = FA1[arr_qpos1] + FA2[arr_qpos2];
			} else if(s == HTS_N) {
				set_base(seq, seq_nt16_str[bam_seqi(seq2, plp[i].qpos)], kh_val(hash, k)->qpos);
				PV1[arr_qpos1] = PV2[arr_qpos2];
				FA1[arr_qpos1] = FA2[arr_qpos2];
			} else if(s2 != HTS_N) {
				++*n_disagreed;
				// Disagreed, both aren't N: N the base, set agrees and p values to 0!
				n_base(seq, kh_val(hash, k)->qpos); // if s2 == HTS_N, do nothing.
				PV1[arr_qpos1] = 0u;
				FA1[arr_qpos1] = 0u;
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
		const int32_t arr_qpos1 = arr_qpos(&plp[i]);
		FA1 = (uint32_t *)array_tag(plp[i].b, "FA");
		PV1 = (uint32_t *)array_tag(plp[i].b, "PV");
		if(FA1[arr_qpos1] < aux->minFA || PV1[arr_qpos1] < aux->minPV ||
			(float)FA1[arr_qpos1] / bam_aux2i(bam_aux_get(plp[i].b, "FM")) < aux->minFR) {
			LOG_DEBUG("Note: PV1[arr_qpos1] value (%u) is greater than minPV now. (%u)\n", PV1[arr_qpos1], aux->minPV);
			continue;
		}
		seq = bam_get_seq(plp[i].b);
		if(bam_seqi(seq, plp[i].qpos) == allele) { // Match!
			++count;
			if((drdata = bam_aux_get(plp[i].b, "DR")) != NULL && bam_aux2i(drdata)) {
				++duplex; // Has DR tag and its value is nonzero.
			}
			if((tmptag = bam_aux_get(plp[i].b, "KR")) != NULL) {
				++overlap;
				bam_aux_del(plp[i].b, tmptag);
			}
		}
	}
	kh_destroy(names, hash);
	*n_duplex = duplex;
	*n_uniobs = count;
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
	int ret = bcf_itr_next(aux->vcf_fp, vcf_iter, vrec);
	return ret;
}

std::vector<std::pair<khint_t, khiter_t>> make_sorted_keys(khash_t(bed) *h) {
	std::vector<std::pair<khint_t, khiter_t>> keyset;
	for(khiter_t ki = kh_begin(aux->bed); ki != kh_end(h); ++ki) {
		if(kh_exist(h, ki)) keyset.push_back(std::pair<khint_t, khiter_t>(kh_key(h, ki), ki));
	}
	__gnu_parallel::sort(keyset.begin(), keyset.end(), [](std::pair<khint_t, khiter_t> p1, std::pair<khint_t, khiter_t> p2) {
		return p1.first < p2.first;
	});
	return keyset;
}

int vet_core(aux_t *aux) {
	int n_plp, vcf_iter_ret;
	const bam_pileup1_t *plp;
	tbx_t *vcf_idx = NULL;
	hts_idx_t *bcf_idx = NULL;
	hts_idx_t *idx = sam_index_load(aux->fp, aux->fp->fn);
	switch(hts_get_format(aux->vcf_fp)->format) {
	case vcf:
		vcf_idx = tbx_index_load(aux->vcf_fp->fn);
		if(!vcf_idx) LOG_ERROR("Could not load TBI index for %s. Indexed vcf required!\n", aux->vcf_fp->fn);
		break;
	case bcf:
		bcf_idx = bcf_index_load(aux->vcf_fp->fn);
		if(!bcf_idx) LOG_ERROR("Could not load CSI index: %s\n", aux->vcf_fp->fn);
		break;
	default:
		LOG_ERROR("Unrecognized variant file type! (%i).\n", hts_get_format(aux->vcf_fp)->format);
		break; // This never happens -- LOG_ERROR exits.
	}
	/*
	if(!(vcf_idx || bcf_idx)) {
		LOG_ERROR("Require an indexed variant file. Abort!\n");
	}
	*/
	bcf1_t *vrec = bcf_init();
	// Unpack all shared data -- up through INFO, but not including FORMAT
	vrec->max_unpack = BCF_UN_FMT;
	vrec->rid = -1;
	hts_itr_t *vcf_iter = NULL;

	std::vector<int32_t> pass_values(DEFAULT_MAX_ALLELES, 0);
	std::vector<int32_t> uniobs_values(DEFAULT_MAX_ALLELES, 0);
	std::vector<int32_t> duplex_values(DEFAULT_MAX_ALLELES, 0);
	std::vector<int32_t> overlap_values(DEFAULT_MAX_ALLELES, 0);
	std::vector<int32_t> fail_values(DEFAULT_MAX_ALLELES, 0);
	bam_plp_t pileup = bam_plp_init(read_bam, (void *)aux);
	bam_plp_set_maxcnt(pileup, max_depth);

	std::vector<std::pair<khint_t, khiter_t>> keys = make_sorted_keys(aux->bed);
	for(auto key: keys) {
		khiter_t ki = key.second;
		for(unsigned j = 0; j < kh_val(aux->bed, ki).n; ++j) {
			int tid, start, stop, pos = -1;

			// Handle coordinates
			tid = kh_key(aux->bed, ki);
			// rid is set to -1 before use. This won't be triggered.
			start = get_start(kh_val(aux->bed, ki).intervals[j]);
			stop = get_stop(kh_val(aux->bed, ki).intervals[j]);
			LOG_DEBUG("Beginning to work through region #%i on contig %s:%i-%i.\n", j + 1, aux->header->target_name[tid], start, stop);

			// Fill vcf_iter from tbi or csi index. If both are null, go through the full file.
			vcf_iter = vcf_idx ? hts_itr_query(vcf_idx->idx, tid, start, stop, tbx_readrec): bcf_itr_queryi(bcf_idx, tid, start, stop);
			LOG_DEBUG("vcf_iter: %p. vcf_idx: %p. bcf_idx: %p.\n", (void *)vcf_iter, (void *)vcf_idx, (void *)bcf_idx);
			if(!vcf_iter) {
				LOG_WARNING("Could not access vcf index. Continuing to next region.\n");
				continue;
			}

			int n_disagreed = 0;
			int n_overlapped = 0;
			int n_duplex = 0;
			while((vcf_iter_ret = bcf_itr_next(aux->vcf_fp, vcf_iter, vrec)) >= 0) {
				if(!bcf_is_snp(vrec)) {
					LOG_DEBUG("Variant isn't a snp. Skip!\n");
					bcf_write(aux->vcf_ofp, aux->vcf_header, vrec);
					continue; // Only handle simple SNVs
				}
				if(!vcf_bed_test(vrec, aux->bed)) {
					LOG_DEBUG("Outside of bed region.\n");
					bcf_write(aux->vcf_ofp, aux->vcf_header, vrec);
					continue; // Only handle simple SNVs
				}
				LOG_DEBUG("Hey, I'm evaluating a variant record now.\n");
				bcf_unpack(vrec, BCF_UN_STR); // Unpack the allele fields
				if (aux->iter) hts_itr_destroy(aux->iter);
				aux->iter = sam_itr_queryi(idx, vrec->rid, vrec->pos, stop);
				plp = bam_plp_auto(pileup, &tid, &pos, &n_plp);
				while ((tid < vrec->rid || (pos < vrec->pos && tid == vrec->rid)) && ((plp = bam_plp_auto(pileup, &tid, &pos, &n_plp)) != 0)) {
					/* Zoom ahead until you're at the correct position */
				}
				if(!plp) {
					if(n_plp == -1) {
						LOG_ERROR("Could not make pileup for region %s:%i-%i. n_plp: %i, pos%i, tid%i.\n", aux->header->target_name[tid], start, stop, n_plp, pos, tid);
					}
					else if(n_plp == 0){
						LOG_WARNING("No reads at position. Skip this variant?\n");
					}
				}
				LOG_DEBUG("tid: %i. rid: %i. var pos: %i.\n", tid, vrec->rid, vrec->pos);
				if(pos != vrec->pos || tid != vrec->rid) {
					LOG_DEBUG("BAM: pos: %i. Contig: %s.\n", pos, aux->header->target_name[tid]);
					LOG_WARNING("Position %s:%i (1-based) not found in pileups in bam. Writing unmodified. Super weird...\n", aux->header->target_name[vrec->rid], vrec->pos + 1);
					bcf_write(aux->vcf_ofp, aux->vcf_header, vrec);
					continue;
				}
				// Check each variant

				// Set i to 1 to skip reference base verification.
				/*bmf_var_tests(bcf1_t *vrec, const bam_pileup1_t *plp, int n_plp, aux_t *aux, std::vector<int>& pass+values,
						std::vector<int>& n_obs, std::vector<int>& n_duplex, std::vector<int>& n_overlaps, std::vector<int> &n_failed,
						int *n_all_overlaps, int *n_all_duplex, int *n_all_disagreed) std::vector<int>& n_obs, std::vector<int>& n_duplex, std::vector<int>& n_overlaps, std::vector<int> &n_failed,
		std::vector<int>& n_pass,*/
				// Reset vectors for each pass.
				std::fill(uniobs_values.begin(), uniobs_values.end(), 0);
				std::fill(duplex_values.begin(), duplex_values.end(), 0);
				std::fill(overlap_values.begin(), overlap_values.end(), 0);
				bmf_var_tests(vrec, plp, n_plp, aux, pass_values, uniobs_values, duplex_values, overlap_values,
							fail_values, n_overlapped, n_duplex, n_disagreed);
				LOG_DEBUG("Adding disc_overlap.\n");
				bcf_update_info_int32(aux->vcf_header, vrec, "DISC_OVERLAP", (void *)&n_disagreed, 1);
				LOG_DEBUG("Adding n_overlap.\n");
				bcf_update_info_int32(aux->vcf_header, vrec, "OVERLAP", (void *)&n_overlapped, 1);
				bcf_update_info_int32(aux->vcf_header, vrec, "DUPLEX_DEPTH", (void *)&n_duplex, 1);
				LOG_DEBUG("n allele: %i.\n", vrec->n_allele);
				bcf_update_info(aux->vcf_header, vrec, "BMF_VET", (void *)(&pass_values[0]), vrec->n_allele, BCF_HT_INT);
				bcf_update_info(aux->vcf_header, vrec, "BMF_FAIL", (void *)(&fail_values[0]), vrec->n_allele, BCF_HT_INT);
				bcf_update_info(aux->vcf_header, vrec, "BMF_DUPLEX", (void *)&duplex_values[0], vrec->n_allele, BCF_HT_INT);
				bcf_update_info(aux->vcf_header, vrec, "BMF_UNIOBS", (void *)&uniobs_values[0], vrec->n_allele, BCF_HT_INT);

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

	if(optind + 1 >= argc) vetter_error("Insufficient arguments. Input bam required!\n", EXIT_FAILURE);
	// Check for required tags.
	if(aux.minAF) check_bam_tag_exit(argv[optind + 1], "AF");
	check_bam_tag_exit(argv[optind + 1], "FA");
	check_bam_tag_exit(argv[optind + 1], "FM");
	check_bam_tag_exit(argv[optind + 1], "FP");
	check_bam_tag_exit(argv[optind + 1], "PV");
	check_bam_tag_exit(argv[optind + 1], "RV");



	strcpy(vcf_wmode, output_bcf ? "wb": "w");
	if(!outvcf) outvcf = strdup("-");
	if(strcmp(outvcf, "-") == 0) {
		LOG_INFO("Emitting to stdout in %s format.\n", output_bcf ? "bcf": "vcf");
	}
	// Open bam
	aux.fp = sam_open_format(argv[optind + 1], "r", &open_fmt);
	if(!aux.fp) LOG_ERROR("Could not open input bam %s. Abort!\n", argv[optind + 1]);
	aux.header = sam_header_read(aux.fp);

	// Open input vcf
	if(!aux.header || aux.header->n_targets == 0) {
		LOG_ERROR("Could not read header from bam %s. Abort!\n", argv[optind + 1]);
	}
	// Open bed file
	// if no bed provided, do whole genome.
	if(bed) aux.bed = parse_bed_hash(bed, aux.header, padding);
	else LOG_ERROR("No bed file provided. Required. Abort!\n");

	if((aux.vcf_fp = vcf_open(argv[optind], "r")) == NULL) LOG_ERROR("Could not open input vcf (%s).\n", argv[optind]);
	if((aux.vcf_header = bcf_hdr_read(aux.vcf_fp)) == NULL) LOG_ERROR("Could not read variant header from file (%s).\n", aux.vcf_fp->fn);

	// Add lines to header
	for(unsigned i = 0; i < COUNT_OF(bmf_header_lines); ++i) {
		LOG_DEBUG("Adding header line %s.\n", bmf_header_lines[i]);
		if(bcf_hdr_append(aux.vcf_header, bmf_header_lines[i])) LOG_ERROR("Could not add header line '%s'. Abort!\n", bmf_header_lines[i]);
	}
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
	if(!aux.vcf_ofp) LOG_ERROR("Could not open output vcf '%s' for writing. Abort!\n", outvcf);
	bcf_hdr_write(aux.vcf_ofp, aux.vcf_header);

	// Open out vcf
	int ret = vet_core(&aux);
	if(ret) {
		fprintf(stderr, "[E:%s:%d] vet_core returned non-zero exit status '%i'. Abort!\n",
				__func__, __LINE__, ret);
	}
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
