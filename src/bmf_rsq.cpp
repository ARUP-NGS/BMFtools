/*  bam_rsq.c */
#include "bmf_rsq.h"

typedef struct rsq_settings {
	FILE *fqh;
	samFile *in;
	samFile *out;
	int cmpkey; // 0 for pos, 1 for unclipped start position
	int mmlim; // Mismatch failure threshold.
	int read_hd_threshold;
	int is_se;
	bam_hdr_t *hdr; // BAM header
	stack_fn fn;
	std::unordered_map<char *, std::string> realign_pairs;
} rsq_settings_t;

void resize_stack(tmp_stack_t *stack, size_t n) {
	if(n > stack->max) {
		stack->max = n;
		stack->a = (bam1_t **)realloc(stack->a, sizeof(bam1_t *) * n);
		if(!stack->a) LOG_EXIT("Failed to reallocate memory for %lu bam1_t * objects. Abort!\n", stack->max);
	} else if(n < stack->n){
		for(uint64_t i = stack->n;i > n;) free(stack->a[--i]->data);
		stack->max = n;
		stack->a = (bam1_t **)realloc(stack->a, sizeof(bam1_t *) * n);
	}
}

std::string bam2cppstr(bam1_t *b)
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
	ksprintf(&ks, "\tFM:i:%i\tFP:i:%i\tNC:i:%i", bam_aux2i(bam_aux_get(b, (char *)"FM")),
			bam_aux2i(bam_aux_get(b, (char *)"FP")),
			bam_aux2i(bam_aux_get(b, (char *)"NC")));
	if((rvdata = bam_aux_get(b, (char *)"RV")) != NULL)
		ksprintf(&ks, "\tRV:i:%i", bam_aux2i(rvdata));
	kputc('\n', &ks);
	seq = bam_get_seq(b);
	seqbuf = (char *)malloc(b->core.l_qseq + 1);
	for (i = 0; i < b->core.l_qseq; ++i) seqbuf[i] = seq_nt16_str[bam_seqi(seq, i)];
	seqbuf[i] = '\0';
	if (b->core.flag & BAM_FREVERSE) { // reverse complement
		for(i = 0; i < b->core.l_qseq>>1; ++i) {
			t = seqbuf[b->core.l_qseq - i - 1];
			seqbuf[b->core.l_qseq - i - 1] = nuc_cmpl(seqbuf[i]);
			seqbuf[i] = nuc_cmpl(t);
		}
		if(b->core.l_qseq&1) seqbuf[i] = nuc_cmpl(seqbuf[i]);
	}
	assert(strlen(seqbuf) == (uint64_t)b->core.l_qseq);
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
	assert(strlen(seqbuf) == (uint64_t)b->core.l_qseq);
	kputs(seqbuf, &ks); kputc('\n', &ks);
	free(seqbuf);
	std::string ret(ks.s);
	free(ks.s);
	return ret;
}

inline int switch_names(char *n1, char *n2) {
	for(;*n1;++n1, ++n2) {
		if(*n1 != *n2) {
			if(*n1 == 'N')
				return 1;
			if(*n2 == 'N')
				return 0;
			return *n1 < *n2;
		}
	}
	return 0; // This never happens.
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
	char *pname = bam_get_qname(p), *bname = bam_get_qname(b);
	if(switch_names(pname, bname)) {
		memcpy(bam_get_qname(p), bam_get_qname(b), b->core.l_qname);
	}
	pFM += bFM;
	bam_aux_del(p, pdata);
	bam_aux_append(p, "FM", 'i', sizeof(int), (uint8_t *)&pFM);
	if((pdata = bam_aux_get(p, "RV")) != NULL) {
		const int pRV = bam_aux2i(pdata) + bam_aux2i(bam_aux_get(b, "RV"));
		bam_aux_del(p, pdata);
		bam_aux_append(p, "RV", 'i', sizeof(int), (uint8_t *)&pRV);
	}
	// Handle NC (Number Changed) tag
	pdata = bam_aux_get(p, "NC");
	bdata = bam_aux_get(b, "NC");
	n_changed = int_tag_zero(pdata) + int_tag_zero(bdata);
	if(pdata) bam_aux_del(p, pdata);

	uint32_t *const bPV = (uint32_t *)array_tag(b, "PV"); // Length of this should be b->l_qseq
	uint32_t *const pPV = (uint32_t *)array_tag(p, "PV"); // Length of this should be b->l_qseq
	uint32_t *const bFA = (uint32_t *)array_tag(b, "FA"); // Length of this should be b->l_qseq
	uint32_t *const pFA = (uint32_t *)array_tag(p, "FA"); // Length of this should be b->l_qseq
	uint8_t *const bSeq = (uint8_t *)bam_get_seq(b);
	uint8_t *const pSeq = (uint8_t *)bam_get_seq(p);
	uint8_t *const bQual = (uint8_t *)bam_get_qual(b);
	uint8_t *const pQual = (uint8_t *)bam_get_qual(p);
	const int qlen = p->core.l_qseq;

	if(p->core.flag & (BAM_FREVERSE)) {
		int qleni1;
		int8_t ps, bs;
		for(int i = 0; i < qlen; ++i) {
			qleni1 = qlen - i - 1;
			ps = bam_seqi(pSeq, qleni1);
			bs = bam_seqi(bSeq, qleni1);
			if(ps == bs) {
				pPV[i] = agreed_pvalues(pPV[i], bPV[i]);
				pFA[i] += bFA[i];
				if(bQual[qleni1] > pQual[qleni1]) pQual[qleni1] = bQual[qleni1];
			} else if(ps == HTS_N) {
				bam_set_base(pSeq, bSeq, qleni1);
				pFA[i] = bFA[i];
				pPV[i] = bPV[i];
				pQual[qleni1] = bQual[qleni1];
				++n_changed; // Note: goes from N to a useable nucleotide.
				continue;
			} else if(bs == HTS_N) {
				continue;
			} else {
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
		int8_t ps, bs;
		for(int i = 0; i < qlen; ++i) {
			ps = bam_seqi(pSeq, i);
			bs = bam_seqi(bSeq, i);
			if(ps == bs) {
				pPV[i] = agreed_pvalues(pPV[i], bPV[i]);
				pFA[i] += bFA[i];
				if(bQual[i] > pQual[i]) pQual[i] = bQual[i];
			} else if(ps == HTS_N) {
				bam_set_base(pSeq, bSeq, i);
				pFA[i] = bFA[i];
				pPV[i] = bPV[i];
				++n_changed; // Note: goes from N to a useable nucleotide.
				continue;
			} else if(bs == HTS_N) {
				continue;
			} else {
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
			} else if((uint32_t)(pQual[i]) > pPV[i]) pQual[i] = (uint8_t)pPV[i];
		}
	}
	bam_aux_append(p, "NC", 'i', sizeof(int), (uint8_t *)&n_changed);
}


void write_stack(tmp_stack_t *stack, rsq_settings_t *settings)
{
	for(unsigned i = 0; i < stack->n; ++i) {
		if(stack->a[i]) {
			uint8_t *data;
			char *qname;
			if((data = bam_aux_get(stack->a[i], "NC")) != NULL) {
				qname = bam_get_qname(stack->a[i]);
				if(settings->realign_pairs.find(qname) == settings->realign_pairs.end())
					settings->realign_pairs[qname] = std::string(bam2cppstr(stack->a[i]));
				else {
					if(stack->a[i]->core.flag & BAM_FREAD2)
						settings->realign_pairs[qname] += bam2cppstr(stack->a[i]);
					else settings->realign_pairs[qname].insert(0, bam2cppstr(stack->a[i]));
				}
			} else sam_write1(settings->out, settings->hdr, stack->a[i]);
			bam_destroy1(stack->a[i]);
			stack->a[i] = NULL;
		}
	}
}


static inline int hd_linear(bam1_t *a, bam1_t *b, int mmlim)
{
	char *aname = (char *)bam_get_qname(a);
	char *bname = (char *)bam_get_qname(b);
	int hd = 0;
	while(*aname)
		if(*aname++ != *bname++)
			if(++hd > mmlim)
				return 0;
	return hd;
}

static inline void flatten_stack_linear(tmp_stack_t *stack, rsq_settings_t *settings)
{
	std::sort(stack->a, stack->a + stack->n, [](const bam1_t *a, const bam1_t *b) {
		if(a) return b ? strcmp(bam_get_qname(a), bam_get_qname(b)) < 0: 0;
		else return b ? 1: 0;
		// Returns 0 if comparing two nulls, and returns true that a NULL lt a valued name
		// Compares strings otherwise.
	});
	for(unsigned i = 0; i < stack->n; ++i) {
		for(unsigned j = i + 1; j < stack->n; ++j) {
			if(hd_linear(stack->a[i], stack->a[j], settings->mmlim) &&
					read_pass_hd(stack->a[i], stack->a[j], settings->read_hd_threshold)) {
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

static inline void rsq_core(rsq_settings_t *settings, tmp_stack_t *stack)
{
	bam1_t *b = bam_init1();
	if(sam_read1(settings->in, settings->hdr, b) < 0)
		LOG_EXIT("Failed to read first record in bam file. Abort!\n");
	while(b->core.flag & (BAM_FSUPPLEMENTARY | BAM_FSECONDARY)) {
		sam_write1(settings->out, settings->hdr, b);
		if(sam_read1(settings->in, settings->hdr, b))
			LOG_EXIT("Could not read first primary alignment in bam (%s). Abort!\n", settings->in->fn);
	}
	stack_insert(stack, b);
	while (LIKELY(sam_read1(settings->in, settings->hdr, b) >= 0)) {
		if(b->core.flag & (BAM_FSUPPLEMENTARY | BAM_FSECONDARY)) {
			sam_write1(settings->out, settings->hdr, b);
			continue;
		}
		if(settings->fn(b, *stack->a) == 0) {
			flatten_stack_linear(stack, settings); // Change this later if the chemistry necessitates it.
			write_stack(stack, settings);
			stack->n = 1;
			stack->a[0] = bam_dup1(b);
		} else stack_insert(stack, b);
	}
	flatten_stack_linear(stack, settings); // Change this later if the chemistry necessitates it.
	write_stack(stack, settings);
	stack->n = 1;
	bam_destroy1(b);
	for(auto& pair: settings->realign_pairs) {
		fprintf(settings->fqh, pair.second.c_str());
	}
}

void bam_rsq_bookends(rsq_settings_t *settings)
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
	if(!stack.a)
		LOG_EXIT("Failed to start array of bam1_t structs...\n");
	rsq_core(settings, &stack); // Core
	free(stack.a);
}

int rsq_usage(int retcode)
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage:  bmftools rsq [-srtu] -f <to_realign.fq> <input.srt.bam> <output.bam>\n\n");
	fprintf(stderr, "Flags:\n"
					"-t	  Mismatch limit. Default: 2\n");
	fprintf(stderr, "-l	  Set bam compression level. Valid: 0-9. (0 == uncompresed)\n");
	fprintf(stderr, "-u	  Flag to use unclipped start positions instead of pos/mpos for identifying potential duplicates.\n"
					"Note: This requires pre-processing with bmftools mark_unclipped.\n");
	return retcode;
}


int rsq_main(int argc, char *argv[])
{
	int c;
	char wmode[4] = {'w', 'b', 0, 0};

	rsq_settings_t settings = {0};
	settings.mmlim = 2;
	settings.cmpkey = POSITION;
	settings.read_hd_threshold = -1;

	char *fqname = NULL;

	if(argc < 3) rsq_usage(EXIT_FAILURE);

	while ((c = getopt(argc, argv, "l:f:t:au?h")) >= 0) {
		switch (c) {
		case 'u': settings.cmpkey = UNCLIPPED; break;
		case 't': settings.mmlim = atoi(optarg); break;
		case 'f': fqname = optarg; break;
		case 'l': wmode[2] = atoi(optarg) + '0'; if(wmode[2] > '9') wmode[2] = '9'; break;
		case 'h': /* fall-through */
		case '?': return rsq_usage(EXIT_SUCCESS);
		}
	}
	if (optind + 2 > argc)
		return rsq_usage(EXIT_FAILURE);
	if(settings.read_hd_threshold < 0) {
		settings.read_hd_threshold = READ_HD_LIMIT;
		LOG_INFO("Unset read HD threshold. Setting to default (%i).\n", settings.read_hd_threshold);
	}

	if(!fqname) {
		fprintf(stderr, "Fastq path for rescued reads required. Abort!\n");
		return rsq_usage(EXIT_FAILURE);
	}

	settings.fqh = fopen(fqname, "w");

	if(!settings.fqh)
		LOG_EXIT("Failed to open output fastq for writing. Abort!\n");


	settings.in = sam_open(argv[optind], "r");
	settings.hdr = sam_hdr_read(settings.in);

	if (settings.hdr == NULL || settings.hdr->n_targets == 0)
		LOG_EXIT("input SAM does not have header. Abort!\n");

	settings.out = sam_open(argv[optind+1], wmode);
	if (settings.in == 0 || settings.out == 0)
		LOG_EXIT("fail to read/write input files\n");
	sam_hdr_write(settings.out, settings.hdr);

	bam_rsq_bookends(&settings);
	bam_hdr_destroy(settings.hdr);
	sam_close(settings.in); sam_close(settings.out);
	if(settings.fqh) fclose(settings.fqh);
	LOG_INFO("Successfully completed bmftools rsq.\n");
	return EXIT_SUCCESS;
}
