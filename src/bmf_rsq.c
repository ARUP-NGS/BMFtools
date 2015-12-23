/*  bam_rsq.c -- duplicate read detection.

	Copyright (C) 2009, 2015 Genome Research Ltd.
	Portions copyright (C) 2009 Broad Institute.

	Author: Heng Li <lh3@sanger.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */
#include "bmf_rsq.h"

void resize_stack(tmp_stack_t *stack, size_t n) {
	if(n > stack->max) {
		stack->max = n;
		stack->a = (bam1_t **)realloc(stack->a, sizeof(bam1_t *) * n);
#if !NDEBUG
			if(!stack->a) {
				fprintf(stderr, "[E:%s] Failed to reallocate memory for %i bam1_t * objects. Abort!\n", __func__, stack->max);
				exit(EXIT_FAILURE);
			}
#endif
	}
	else if(n < stack->n){
		for(uint64_t i = stack->n;i > n;) bam_destroy1(stack->a[--i]);
		stack->max = n;
		stack->a = (bam1_t **)realloc(stack->a, sizeof(bam1_t *) * n);
	}
}

static inline void update_bam1(bam1_t *p, bam1_t *b)
{
	uint8_t *bdata, *pdata;
	uint32_t *const bPV = (uint32_t *)array_tag(b, "PV"); // Length of this should be b->l_qseq
	uint32_t *const pPV = (uint32_t *)array_tag(p, "PV"); // Length of this should be b->l_qseq
	uint32_t *const bFA = (uint32_t *)array_tag(b, "FA"); // Length of this should be b->l_qseq
	uint32_t *const pFA = (uint32_t *)array_tag(p, "FA"); // Length of this should be b->l_qseq
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

	// Check for required PV and FA tags
	if(!bPV || !pPV) {
		fprintf(stderr, "Required PV tag not found. Abort mission! Read names: %s, %s.\n", bam_get_qname(b), bam_get_qname(p));
		exit(EXIT_FAILURE);
	}
	if(!bFA || !pFA) {
		fprintf(stderr, "Required FA tag not found. Abort mission!\n");
		exit(EXIT_FAILURE);
	}

	uint8_t *const bSeq = (uint8_t *)bam_get_seq(b);
	uint8_t *const pSeq = (uint8_t *)bam_get_seq(p);
	uint8_t *const bQual = (uint8_t *)bam_get_qual(b);
	uint8_t *const pQual = (uint8_t *)bam_get_qual(p);
	if(!(bSeq && pSeq && bQual && pQual)) {
		fprintf(stderr, "Qual strings or sequence strings are null. Abort!\n");
	}
	const int qlen = p->core.l_qseq;
	int qleni1;
	if(p->core.flag & (BAM_FREVERSE)) {
		for(int i = 0; i < qlen; ++i) {
			qleni1 = qlen - i - 1;
			if(bam_seqi(pSeq, qleni1) == bam_seqi(bSeq, qleni1)) {
				pPV[i] = agreed_pvalues(pPV[i], bPV[i]);
				pFA[i] += bFA[i];
				if(bQual[qleni1] > pQual[qleni1])
					pQual[qleni1] = bQual[qleni1];
			}
			else if(bam_seqi(pSeq, qleni1) == HTS_N) {
				set_base(pSeq, bSeq, qleni1);
				pFA[i] = bFA[i];
				pPV[i] = bPV[i];
				pQual[qleni1] = bQual[qleni1];
				++n_changed; // Note: goes from N to a useable nucleotide.
				continue;
			} else {
				pPV[i] = (pPV[i] > bPV[i]) ? disc_pvalues(pPV[i], bPV[i]) : disc_pvalues(bPV[i], pPV[i]);
				if(bam_seqi(bSeq, qleni1) != HTS_N)
					set_base(pSeq, bSeq, qleni1);
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
			if((uint32_t)(pQual[qleni1]) > pPV[i]) pQual[qleni1] = pPV[i];
		}
	}
	else {
		for(int i = 0; i < qlen; ++i) {
			if(bam_seqi(pSeq, i) == bam_seqi(bSeq, i)) {
				pPV[i] = agreed_pvalues(pPV[i], bPV[i]);
				pFA[i] += bFA[i];
				if(bQual[i] > pQual[i])
					pQual[i] = bQual[i];
			}
			else if(bam_seqi(pSeq, i) == HTS_N) {
				set_base(pSeq, bSeq, i);
				pFA[i] = bFA[i];
				pPV[i] = bPV[i];
				++n_changed; // Note: goes from N to a useable nucleotide.
				continue;
			}
			else {
				pPV[i] = (pPV[i] > bPV[i]) ? disc_pvalues(pPV[i], bPV[i]) : disc_pvalues(bPV[i], pPV[i]);
				if(bam_seqi(bSeq, i) != HTS_N)
					set_base(pSeq, bSeq, i);
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
			if((uint32_t)(pQual[i]) > pPV[i]) pQual[i] = pPV[i];
		}
	}
	bam_aux_append(p, "NC", 'i', sizeof(int), (uint8_t *)&n_changed);
}

void bam2ffq(bam1_t *b, FILE *fp)
{
	int i;
	int qlen = b->core.l_qseq;
	char seqbuf[SEQBUF_SIZE];
	uint8_t *seq = bam_get_seq(b);
	for (i = 0; i < qlen; ++i)
		seqbuf[i] = bam_seqi(seq, i);
	if (b->core.flag & BAM_FREVERSE) { // reverse complement
		for (i = 0; i < qlen>>1; ++i) {
			int8_t t = seq_comp_table[(int8_t)seqbuf[qlen - 1 - i]];
			seqbuf[qlen - 1 - i] = seq_comp_table[(int8_t)seqbuf[i]];
			seqbuf[i] = t;
		}
		if (qlen&1) seqbuf[i] = seq_comp_table[(int8_t)seqbuf[i]];
	}
	for (i = 0; i < qlen; ++i)
		seqbuf[i] = seq_nt16_str[(int8_t)seqbuf[i]];
	seqbuf[qlen] = '\0';
#if !NDEBUG
	fprintf(stderr, "seqbuf: %s.\n", seqbuf);
#endif
	char comment[3000] = "";
	uint32_t *pv = (uint32_t *)array_tag(b, (char *)"PV");
	uint32_t *fa = (uint32_t *)array_tag(b, (char *)"FA");
	append_csv_buffer(b->core.l_qseq, pv, comment, (char *)"PV:B:I");
	strcat(comment, "\t");
	append_csv_buffer(b->core.l_qseq, fa, comment, (char *)"FA:B:I");
	append_int_tag(comment, (char *)"FM", bam_aux2i(bam_aux_get(b, (char *)"FM")));
	const uint8_t *rvdata = bam_aux_get(b, (char *)"RV");
	if(rvdata)
		append_int_tag(comment, (char *)"RV", bam_aux2i(rvdata));
	append_int_tag(comment, (char *)"FP", bam_aux2i(bam_aux_get(b, (char *)"FP")));
	append_int_tag(comment, (char *)"NC", bam_aux2i(bam_aux_get(b, (char *)"NC")));
#if !NDEBUG
	fprintf(stderr, "comment string: %s.\n", comment);
#endif
	fprintf(fp, "@%s %s\n%s\n+\n", (char *)bam_get_qname(b), comment, seqbuf);
	char *qual = (char *)bam_get_qual(b);
	for(i = 0; i < qlen; ++i)
		seqbuf[i] = 33 + qual[i];
	if (b->core.flag & BAM_FREVERSE) { // reverse
		for (i = 0; i < qlen>>1; ++i) {
			const int8_t t = seqbuf[qlen - 1 - i];
			seqbuf[qlen - 1 - i] = seqbuf[i];
			seqbuf[i] = t;
		}
	}
	fprintf(fp, "%s\n", seqbuf);
	return;
}



void write_stack(tmp_stack_t *stack, pr_settings_t *settings)
{
	for(int i = 0; i < stack->n; ++i) {
		if(stack->a[i]) {
			uint8_t *data;
			if((data = bam_aux_get(stack->a[i], "NC")) != NULL) {
				if(bam_aux2i(data) == 0)
					sam_write1(settings->out, settings->hdr, stack->a[i]);
				else
#if !NDEBUG
				{ fprintf(stderr, "About to write a record to the fastq!");
#endif
					bam2ffq(stack->a[i], settings->fqh);
#if !NDEBUG
				}
#endif
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


void bam_rsqse_core(pr_settings_t *settings)
{
	/*
	bam1_t *b;
	tmp_stack_t stack;
	resize_stack(&stack, STACK_START);
	memset(&stack, 0, sizeof(tmp_stack_t) * stack->max);

	b = bam_init1();
	while (sam_read1(settings.in, settings.hdr, b) >= 0) {
		if(stack.n) {

		}
	}
	bam_destroy1(b);
	*/
}

static inline int hd_linear(bam1_t *a, bam1_t *b, int mmlim)
{
	char *aname = (char *)bam_get_qname(a);
	char *bname = (char *)bam_get_qname(b);
	int l_qname = a->core.l_qname - 1; // Skip the terminal null in comparison.
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
	for(int i = 0; i < stack->n; ++i) {
		for(int j = i + 1; j < stack->n; ++j) {
			if(hd_linear(stack->a[i], stack->a[j], settings->mmlim)) {
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

static inline void rsq_core_pos(pr_settings_t *settings, tmp_stack_t *stack)
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
		if(same_stack_pos(b, *stack->a)) {
			stack_insert(stack, b);
			continue;
		} else {
			flatten_stack_linear(stack, settings); // Change this later if the chemistry necessitates it.
			write_stack(stack, settings);
			stack->n = 1;
			stack->a[0] = bam_dup1(b);
		}
	}
	bam_destroy1(b);
}

static inline void rsq_core_ucs(pr_settings_t *settings, tmp_stack_t *stack)
{
	bam1_t *b = bam_init1();
	if(sam_read1(settings->in, settings->hdr, b) < 0) {
		fprintf(stderr, "[E:%s] Failed to read first record in bam file. Abort!\n", __func__);
		exit(EXIT_FAILURE);
	}
	while(b->core.flag & (BAM_FSUPPLEMENTARY | BAM_FSECONDARY)) sam_write1(settings->out, settings->hdr, b);
	stack_insert(stack, b);
	if(!(settings->in && settings->hdr)) {
		fprintf(stderr, "[E:%s] Failed to open input bam.\n", __func__);
		exit(EXIT_FAILURE);
	}
	while (LIKELY(sam_read1(settings->in, settings->hdr, b)) >= 0) {
		if(b->core.flag & (BAM_FSUPPLEMENTARY | BAM_FSECONDARY)) {
			sam_write1(settings->out, settings->hdr, b);
			continue;
		}
		if(same_stack_ucs(b, *stack->a)) {
			stack_insert(stack, b);
			continue;
		} else {
			flatten_stack_linear(stack, settings); // Change this later if the chemistry necessitates it.
			write_stack(stack, settings);
			stack->n = 1;
			stack->a[0] = bam_dup1(b);
		}
	}
	bam_destroy1(b);
}


static inline void rsq_core(pr_settings_t *settings, tmp_stack_t *stack)
{
	settings->cmpkey ? rsq_core_ucs(settings, stack): rsq_core_pos(settings, stack);
}


void bam_rsq_bookends(pr_settings_t *settings)
{
	tmp_stack_t stack;
	memset(&stack, 0, sizeof(tmp_stack_t));
	resize_stack(&stack, STACK_START);
	if(!(settings->in && settings->hdr && settings->out)) {
		fprintf(stderr, "Failed to read input/output files....\n");
		exit(EXIT_FAILURE);
	}
	if(!stack.a) {
		fprintf(stderr, "Failed to start array of bam1_t structs...\n");
		exit(EXIT_FAILURE);
	}
	rsq_core(settings, &stack); // Core
	free(stack.a);
}

int pr_usage(void)
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage:  bmftools rsq [-srtu] -f <to_realign.fq> <input.srt.bam> <output.bam>\n\n");
	fprintf(stderr, "Flags:\n"
					"-s      Rescue for SE reads [Not implemented]\n");
	fprintf(stderr, "-r      Realign reads with no changed bases. Default: False.\n");
	fprintf(stderr, "-t      Mismatch limit. Default: 2\n");
	fprintf(stderr, "-l      Set bam compression level. Valid: 0-9. (0 == uncompresed)\n");
	fprintf(stderr, "-u      Flag to use unclipped start positions instead of pos/mpos for identifying potential duplicates.\n"
					"Note: This requires pre-processing with bmftools mark_unclipped.\n");
	return 1;
}


int rsq_main(int argc, char *argv[])
{
	int c;
	char wmode[4] = {'w', 'b', 0, 0};
	sam_global_args ga = SAM_GLOBAL_ARGS_INIT;

	static const struct option lopts[] = {
		SAM_OPT_GLOBAL_OPTIONS('-', 0, 0, 0, 0),
		{ NULL, 0, NULL, 0 }
	};

	pr_settings_t *settings = (pr_settings_t *)malloc(sizeof(pr_settings_t));
	settings->fqh = NULL;
	settings->in = NULL;
	settings->out = NULL;
	settings->hdr = NULL;
	settings->mmlim = 2;
	settings->cmpkey = POSITION;
	settings->is_se = 0;
	settings->realign_unchanged = 0;

	char fqname[200] = "";

	while ((c = getopt_long(argc, argv, "l:f:t:aur?h", lopts, NULL)) >= 0) {
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

	if(strcmp(fqname, "") == 0) {
		fprintf(stderr, "Fastq path for rescued reads required. Abort!\n");
		return pr_usage();
	}

	settings->fqh = fopen(fqname, "w");

	if(!settings->fqh) {
		fprintf(stderr, "Failed to open output fastq for writing. Abort!\n");
		exit(EXIT_FAILURE);
	}


	settings->in = sam_open_format(argv[optind], "rb", &ga.in);
	settings->hdr = sam_hdr_read(settings->in);
	if (settings->hdr == NULL || settings->hdr->n_targets == 0) {
		fprintf(stderr, "[bam_rsq] input SAM does not have header. Abort!\n");
		return 1;
	}

	settings->out = sam_open_format(argv[optind+1], wmode, &ga.out);
	if (settings->in == 0 || settings->out == 0) {
		fprintf(stderr, "[bam_rsq] fail to read/write input files\n");
		return 1;
	}
	sam_hdr_write(settings->out, settings->hdr);

	if (settings->is_se) bam_rsqse_core(settings);
	else bam_rsq_bookends(settings);
	bam_hdr_destroy(settings->hdr);
	sam_close(settings->in); sam_close(settings->out);
	if(settings->fqh)
		fclose(settings->fqh);
	return 0;
}
