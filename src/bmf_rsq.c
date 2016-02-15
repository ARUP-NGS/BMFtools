/*  bam_rsq.c */
#include "bmf_rsq.h"

void resize_stack(tmp_stack_t *stack, size_t n) {
	if(n > stack->max) {
		stack->max = n;
		stack->a = (bam1_t **)realloc(stack->a, sizeof(bam1_t *) * n);
		if(!stack->a) LOG_ERROR("Failed to reallocate memory for %lu bam1_t * objects. Abort!\n", stack->max);
	} else if(n < stack->n){
		for(uint64_t i = stack->n;i > n;) bam_destroy1(stack->a[--i]);
		stack->max = n;
		stack->a = (bam1_t **)realloc(stack->a, sizeof(bam1_t *) * n);
	}
}

void bam2ffq(bam1_t *b, FILE *fp)
{
	char *qual, *seqbuf;
	int i;
	uint8_t *seq, *rvdata;
	uint32_t *pv, *fa;
	int8_t t;
	kstring_t ks = {0, 0, NULL};
	ksprintf(&ks, "@%s PV:B:I", bam_get_qname(b));
	pv = (uint32_t *)array_tag(b, (char *)"PV");
	fa = (uint32_t *)array_tag(b, (char *)"FA");
	for(i = 0; i < b->core.l_qseq; ++i) ksprintf(&ks, ",%u", pv[i]);
	kputs("\tFA:B:I", &ks);
	for(i = 0; i < b->core.l_qseq; ++i) ksprintf(&ks, ",%u", fa[i]);
	seq = bam_get_seq(b);
	for (i = 0; i < b->core.l_qseq; ++i) kputc(bam_seqi(seq, i), &ks);
	ksprintf(&ks, "\tFM:i:%i\tFP:i:%i\tNC:i:%i", bam_aux2i(bam_aux_get(b, (char *)"FM")),
			bam_aux2i(bam_aux_get(b, (char *)"FP")),
			bam_aux2i(bam_aux_get(b, (char *)"NC")));
	if((rvdata = bam_aux_get(b, (char *)"RV")) != NULL)
		ksprintf(&ks, "\tRV:i:%i", bam_aux2i(rvdata));
	kputc('\n', &ks);
	seqbuf = (char *)malloc(b->core.l_qseq + 1);
	if (b->core.flag & BAM_FREVERSE) { // reverse complement
		for(i = 0; i < b->core.l_qseq>>1; ++i) {
			t = seq_comp_table[(int8_t)seqbuf[b->core.l_qseq - 1 - i]];
			seqbuf[b->core.l_qseq - 1 - i] = seq_comp_table[(int8_t)seqbuf[i]];
			seqbuf[i] = t;
		}
		if(b->core.l_qseq&1) seqbuf[i] = seq_comp_table[(int8_t)seqbuf[i]];
	}
	for (i = 0; i < b->core.l_qseq; ++i) seqbuf[i] = seq_nt16_str[(int8_t)seqbuf[i]];
	seqbuf[i] = '\0'; // i == b->core.l_qseq
	kputs(seqbuf, &ks);
	kputs("\n+\n", &ks);
	qual = (char *)bam_get_qual(b);
	for(i = 0; i < b->core.l_qseq; ++i) seqbuf[i] = 33 + qual[i];
	if (b->core.flag & BAM_FREVERSE) { // reverse
		for (i = 0; i < b->core.l_qseq>>1; ++i) {
			t = seqbuf[b->core.l_qseq - 1 - i];
			seqbuf[b->core.l_qseq - 1 - i] = seqbuf[i];
			seqbuf[i] = t;
		}
	}
	kputs(seqbuf, &ks); kputc('\n', &ks);
	free(seqbuf);
	fputs(ks.s, fp);
	free(ks.s);
	return;
}

static inline void update_bam1(bam1_t *p, bam1_t *b)
{
	uint8_t *bdata, *pdata;
	int n_changed;
/*
	if(!b || !p) {
		// If the
		fprintf(stderr, "One of these records is null. Abort!\n");
		exit(EXIT_FAILURE);
	}
*/
	bdata = bam_aux_get(b, "FM");
	pdata = bam_aux_get(p, "FM");
	if(!bdata || !pdata) {
		fprintf(stderr, "Required FM tag not found. Abort mission!\n");
		exit(EXIT_FAILURE);
	}
	int bFM = bam_aux2i(bdata);
	int pFM = bam_aux2i(pdata);
	if(pFM < bFM)
		memcpy(bam_get_qname(p), bam_get_qname(b), b->core.l_qname);
	pFM += bFM;
	bam_aux_del(p, pdata);
	bam_aux_append(p, "FM", 'i', sizeof(int), (uint8_t *)&pFM);
#if !NDEBUG
	assert(pFM == bam_aux2i(bam_aux_get(p, "FM")));
#endif
	pdata = bam_aux_get(p, "RV");
	if(pdata) {
		const int pRV = bam_aux2i(pdata) + bam_aux2i(bam_aux_get(b, "RV"));
		bam_aux_del(p, pdata);
		bam_aux_append(p, "RV", 'i', sizeof(int), (uint8_t *)&pRV);
	}
	// Handle NC (Number Changed) tag
	pdata = bam_aux_get(p, "NC");
	bdata = bam_aux_get(b, "NC");
	if(pdata) {
		n_changed = (bdata) ? bam_aux2i(bdata) + bam_aux2i(pdata): bam_aux2i(pdata);
		bam_aux_del(p, pdata);
	} else n_changed = (bdata) ? bam_aux2i(bdata): 0;

	uint32_t *const bPV = (uint32_t *)array_tag(b, "PV"); // Length of this should be b->l_qseq
	uint32_t *const pPV = (uint32_t *)array_tag(p, "PV"); // Length of this should be b->l_qseq
	uint32_t *const bFA = (uint32_t *)array_tag(b, "FA"); // Length of this should be b->l_qseq
	uint32_t *const pFA = (uint32_t *)array_tag(p, "FA"); // Length of this should be b->l_qseq
#if !NDEBUG
	// Check for required PV and FA tags
	if(!bPV || !pPV) {
		LOG_ERROR("Required PV tag not found. Abort mission! Read names: %s, %s.\n",
				bam_get_qname(b), bam_get_qname(p));
	}
	if(!bFA || !pFA) {
		LOG_ERROR("Required FA tag not found. Abort mission! Read names: %s, %s.\n",
				bam_get_qname(b), bam_get_qname(p));
	}
#endif

	uint8_t *const bSeq = (uint8_t *)bam_get_seq(b);
	uint8_t *const pSeq = (uint8_t *)bam_get_seq(p);
	uint8_t *const bQual = (uint8_t *)bam_get_qual(b);
	uint8_t *const pQual = (uint8_t *)bam_get_qual(p);
#if !NDEBUG
	if(!(bSeq && pSeq && bQual && pQual)) {
		LOG_ERROR("Qual strings or sequence strings are null. Abort!\n");
	}
#endif
	const int qlen = p->core.l_qseq;
	if(p->core.flag & (BAM_FREVERSE)) {
		int qleni1;
		for(int i = 0; i < qlen; ++i) {
			qleni1 = qlen - i - 1;
			if(bam_seqi(pSeq, qleni1) == bam_seqi(bSeq, qleni1)) {
				pPV[i] = agreed_pvalues(pPV[i], bPV[i]);
				pFA[i] += bFA[i];
				if(bQual[qleni1] > pQual[qleni1]) pQual[qleni1] = bQual[qleni1];
			} else if(bam_seqi(pSeq, qleni1) == HTS_N) {
				bam_set_base(pSeq, bSeq, qleni1);
				pFA[i] = bFA[i];
				pPV[i] = bPV[i];
				pQual[qleni1] = bQual[qleni1];
				++n_changed; // Note: goes from N to a useable nucleotide.
				continue;
			} else if(bam_seqi(bSeq, qleni1) == HTS_N) continue;
			else {
				if(pPV[i] > bPV[i]) {
					bam_set_base(pSeq, bSeq, qleni1);
					pPV[i] = disc_pvalues(pPV[i], bPV[i]);
				} else pPV[i] = disc_pvalues(bPV[i], pPV[i]);
				pFA[i] = bFA[i];
				pQual[qleni1] = bQual[qleni1];
				++n_changed;
			}
			if(pPV[i] < 3) {
				pFA[i] = 0;
				pPV[i] = 0;
				pQual[qleni1] = 2;
				n_base(pSeq, qleni1);
				continue;
			}
			if((uint32_t)(pQual[qleni1]) > pPV[i]) pQual[qleni1] = (uint8_t)pPV[i];
		}
	} else {
		for(int i = 0; i < qlen; ++i) {
			if(bam_seqi(pSeq, i) == bam_seqi(bSeq, i)) {
				pPV[i] = agreed_pvalues(pPV[i], bPV[i]);
				pFA[i] += bFA[i];
				if(bQual[i] > pQual[i]) pQual[i] = bQual[i];
			} else if(bam_seqi(pSeq, i) == HTS_N) {
				bam_set_base(pSeq, bSeq, i);
				pFA[i] = bFA[i];
				pPV[i] = bPV[i];
				++n_changed; // Note: goes from N to a useable nucleotide.
				continue;
			} else if(bam_seqi(bSeq, i) == HTS_N) continue;
			else {
				if(pPV[i] > bPV[i]) {
					bam_set_base(pSeq, bSeq, i);
					pPV[i] = disc_pvalues(pPV[i], bPV[i]);
				} else pPV[i] = disc_pvalues(bPV[i], pPV[i]);
				pFA[i] = bFA[i];
				pQual[i] = bQual[i];
				++n_changed;
			}
			if(pPV[i] < 3) {
				pFA[i] = 0;
				pPV[i] = 0;
				pQual[i] = 2;
				n_base(pSeq, i);
				continue;
			}
			if((uint32_t)(pQual[i]) > pPV[i]) pQual[i] = (uint8_t)pPV[i];
		}
	}
	bam_aux_append(p, "NC", 'i', sizeof(int), (uint8_t *)&n_changed);
}


void write_stack(tmp_stack_t *stack, pr_settings_t *settings)
{
	for(unsigned i = 0; i < stack->n; ++i) {
		if(stack->a[i]) {
			uint8_t *data;
			if((data = bam_aux_get(stack->a[i], "NC")) != NULL) {
				if(bam_aux2i(data) == 0)
					sam_write1(settings->out, settings->hdr, stack->a[i]);
				else
					bam2ffq(stack->a[i], settings->fqh);
			}
			else {
#if !NDEBUG
				if(bam_aux_get(stack->a[i], "NC")) {
					fprintf(stderr, "NC: %i.\n", bam_aux2i(bam_aux_get(stack->a[i], "NC")));
				}
#endif
				sam_write1(settings->out, settings->hdr, stack->a[i]);
			}
			bam_destroy1(stack->a[i]);
			stack->a[i] = NULL;
		}
	}
}


static inline int hd_linear(bam1_t *a, bam1_t *b, int mmlim)
{
	char *aname = (char *)bam_get_qname(a);
	char *bname = (char *)bam_get_qname(b);
	unsigned l_qname = a->core.l_qname - 1; // Skip the terminal null in comparison.
	int hd = 0;
	for(uint64_t i = 0; i < l_qname; ++i) {
		if(aname[i] != bname[i]) {
			if(++hd > mmlim)
				return 0;
		}
	}
	return hd;
}

static inline void flatten_stack_linear(tmp_stack_t *stack, pr_settings_t *settings)
{
	for(unsigned i = 0; i < stack->n; ++i) {
		for(unsigned j = i + 1; j < stack->n; ++j) {
			if(hd_linear(stack->a[i], stack->a[j], settings->mmlim) && read_pass_hd(stack->a[i], stack->a[j], settings->read_hd_threshold)) {
				LOG_DEBUG("Flattening record with qname %s into record with qname %s.\n", bam_get_qname(stack->a[i]), bam_get_qname(stack->a[j]));
				update_bam1(stack->a[j], stack->a[i]);
				bam_destroy1(stack->a[i]);
				stack->a[i] = NULL;
				break;
				// "break" in case there are multiple within hamming distance.
				// Otherwise, I'll end up having memory mistakes.
				// Besides, that read set will get merged into the later read in the set.
			}
		}
	}
}

static inline void rsq_core(pr_settings_t *settings, tmp_stack_t *stack)
{
	bam1_t *b = bam_init1();
	if(sam_read1(settings->in, settings->hdr, b) < 0) {
		fprintf(stderr, "[E:%s] Failed to read first record in bam file. Abort!\n", __func__);
		exit(EXIT_FAILURE);
	}
	while(b->core.flag & (BAM_FSUPPLEMENTARY | BAM_FSECONDARY)) sam_write1(settings->out, settings->hdr, b);
	stack_insert(stack, b);
	if(!(settings->in && settings->hdr)) {
		fprintf(stderr, "Failed to open input bam... WTF?\n");
		exit(EXIT_FAILURE);
	}
	while (LIKELY(sam_read1(settings->in, settings->hdr, b) >= 0)) {
		if(b->core.flag & (BAM_FSUPPLEMENTARY | BAM_FSECONDARY)) {
			sam_write1(settings->out, settings->hdr, b);
			continue;
		}
		if(settings->fn(b, *stack->a)) {
			stack_insert(stack, b);
			continue;
		} else {
			flatten_stack_linear(stack, settings); // Change this later if the chemistry necessitates it.
			write_stack(stack, settings);
			stack->n = 1;
			stack->a[0] = bam_dup1(b);
		}
	}
	flatten_stack_linear(stack, settings); // Change this later if the chemistry necessitates it.
	write_stack(stack, settings);
	stack->n = 1;
	bam_destroy1(b);
}

void bam_rsq_bookends(pr_settings_t *settings)
{
	if(settings->is_se) {
		if(settings->cmpkey)
			settings->fn = &same_stack_ucs_se;
		else
			settings->fn = &same_stack_pos_se;
	} else {
		if(settings->cmpkey)
			settings->fn = &same_stack_ucs;
		else
			settings->fn = &same_stack_pos;
	}


	tmp_stack_t stack;
	memset(&stack, 0, sizeof(tmp_stack_t));
	resize_stack(&stack, STACK_START);
	if(!stack.a) {
		fprintf(stderr, "[E:%s] Failed to start array of bam1_t structs...\n", __func__);
		exit(EXIT_FAILURE);
	}
	else rsq_core(settings, &stack); // Core
	free(stack.a);
}

int pr_usage(void)
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage:  bmftools rsq [-srtu] -f <to_realign.fq> <input.srt.bam> <output.bam>\n\n");
	fprintf(stderr, "Flags:\n"
					"-s	  Rescue for SE reads [Not implemented]\n");
	fprintf(stderr, "-r	  Realign reads with no changed bases. Default: False.\n");
	fprintf(stderr, "-t	  Mismatch limit. Default: 2\n");
	fprintf(stderr, "-l	  Set bam compression level. Valid: 0-9. (0 == uncompresed)\n");
	fprintf(stderr, "-u	  Flag to use unclipped start positions instead of pos/mpos for identifying potential duplicates.\n"
					"Note: This requires pre-processing with bmftools mark_unclipped.\n");
	return 1;
}


int rsq_main(int argc, char *argv[])
{
	int c;
	char wmode[4] = {'w', 'b', 0, 0};

	pr_settings_t *settings = (pr_settings_t *)malloc(sizeof(pr_settings_t));
	settings->fqh = NULL;
	settings->in = NULL;
	settings->out = NULL;
	settings->hdr = NULL;
	settings->mmlim = 2;
	settings->cmpkey = POSITION;
	settings->is_se = 0;
	settings->realign_unchanged = 0;
	settings->read_hd_threshold = -1;

	char fqname[200] = "";

	while ((c = getopt(argc, argv, "l:f:t:aur?h")) >= 0) {
		switch (c) {
		case 'r': settings->realign_unchanged = 1; break;
		case 'u': settings->cmpkey = UNCLIPPED; break;
		case 't': settings->mmlim = atoi(optarg); break;
		case 'f': strcpy(fqname, optarg); break;
		case 'l': wmode[2] = atoi(optarg) + '0'; if(wmode[2] > '9') wmode[2] = '9'; break;
		case 'h': /* fall-through */
		case '?': return pr_usage();
		}
	}
	if (optind + 2 > argc)
		return pr_usage();
	if(settings->read_hd_threshold < 0) {
		settings->read_hd_threshold = READ_HD_LIMIT;
		LOG_INFO("Unset read HD threshold. Setting to default (%i).\n", settings->read_hd_threshold);
	}

	if(strcmp(fqname, "") == 0) {
		fprintf(stderr, "Fastq path for rescued reads required. Abort!\n");
		return pr_usage();
	}

	settings->fqh = fopen(fqname, "w");

	if(!settings->fqh) {
		fprintf(stderr, "Failed to open output fastq for writing. Abort!\n");
		exit(EXIT_FAILURE);
	}


	settings->in = sam_open(argv[optind], "r");
	settings->hdr = sam_hdr_read(settings->in);
	if (settings->hdr == NULL || settings->hdr->n_targets == 0) {
		LOG_ERROR("input SAM does not have header. Abort!\n");
	}

	settings->out = sam_open(argv[optind+1], wmode);
	if (settings->in == 0 || settings->out == 0) {
		LOG_ERROR("fail to read/write input files\n");
	}
	sam_hdr_write(settings->out, settings->hdr);

	if(!(settings->in && settings->hdr && settings->out)) {
		LOG_ERROR("Failed to read input/output files....\n");
	}
	bam_rsq_bookends(settings);
	bam_hdr_destroy(settings->hdr);
	sam_close(settings->in); sam_close(settings->out);
	if(settings->fqh)
		fclose(settings->fqh);
	return EXIT_SUCCESS;
}
