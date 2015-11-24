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
#include "bam_rsq.h"
static inline void update_bam1(bam1_t *p, bam1_t *b)
{
	uint8_t *bdata, *pdata;
	uint32_t *bPV = (uint32_t *)array_tag(b, "PV"); // Length of this should be b->l_qseq
	uint32_t *pPV = (uint32_t *)array_tag(p, "PV"); // Length of this should be b->l_qseq
	uint32_t *bFA = (uint32_t *)array_tag(b, "FA"); // Length of this should be b->l_qseq
	uint32_t *pFA = (uint32_t *)array_tag(p, "FA"); // Length of this should be b->l_qseq
	if(!b || !p) {
		// If the
		fprintf(stderr, "One of these records is null. Abort!\n");
		exit(EXIT_FAILURE);
	}
	bdata = bam_aux_get(b, "FM");
	pdata = bam_aux_get(p, "FM");
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
		int pRV = bam_aux2i(pdata);
		pRV += bam_aux2i(bam_aux_get(b, "RV"));
		bam_aux_del(p, pdata);
		bam_aux_append(p, "RV", 'i', sizeof(int), (uint8_t *)&pRV);
	}
	/*
	pdata = bam_aux_get(p, "NN");
	bdata = bam_aux_get(b, "NN");
	int mask;
	if(pdata) {
		if(bdata)
			mask = bam_aux2i(bdata) + bam_aux2i(pdata);
		else
			mask = bam_aux2i(pdata);
		bam_aux_del(p, pdata);
	}
	else if(bdata)
		mask = bam_aux2i(bdata);
	else
		mask = 0;
	*/
	pdata = bam_aux_get(p, "NC");
	bdata = bam_aux_get(b, "NC");
	int n_changed;
	if(pdata) {
		if(bdata)
			n_changed = bam_aux2i(bdata) + bam_aux2i(pdata);
		else
			n_changed = bam_aux2i(pdata);
		bam_aux_del(p, pdata);
	}
	else if(bdata)
		n_changed = bam_aux2i(bdata);
	else
		n_changed = 0;
/*
#if !NDEBUG
	for(int i = 0; i < b->core.l_qseq; ++i) {
		fprintf(stderr, "%"PRIu32",", bPV[i]);
	}
	fprintf(stderr, "\n");
	exit(0);
#endif
*/
	if(!bPV || !pPV) {
		fprintf(stderr, "Required PV tag not found. Abort mission! Read names: %s, %s.\n", bam_get_qname(b), bam_get_qname(p));
		exit(EXIT_FAILURE);
	}
	if(!bFA || !pFA) {
		fprintf(stderr, "Required FA tag not found. Abort mission!\n");
		exit(EXIT_FAILURE);
	}

	uint8_t *bSeq = (uint8_t *)bam_get_seq(b);
	uint8_t *pSeq = (uint8_t *)bam_get_seq(p);
	uint8_t *bQual = (uint8_t *)bam_get_qual(b);
	uint8_t *pQual = (uint8_t *)bam_get_qual(p);
	if(!(bSeq && pSeq && bQual && pQual)) {
		fprintf(stderr, "Qual strings or sequence strings are null. Abort!\n");
	}
	int qlen = p->core.l_qseq;
	if(p->core.flag & (BAM_FREVERSE)) {
		for(int i = 0; i < qlen; ++i) {
			const int qleni1 = qlen - i - 1;
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
			}
			else {
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

static inline int annealed_check(bam1_t *a, bam1_t *b, int mmlim)
{
	char *aname = (char *)bam_get_qname(a);
	char *bname = (char *)bam_get_qname(b);
	int lq_2 = (a->core.l_qname - 1) / 2;
	int hd = 0;
	int i;
	for(i = 0; i < lq_2; ++i) {
		if(aname[i] != bname[i + lq_2]) {
			if(++hd > mmlim) {
				return 0;
			}
		}
		if(aname[i + lq_2] != bname[i]) {
			if(++hd > mmlim) {
				return 0;
			}
		}
	}
	return hd * -1;
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

static inline void pr_loop_pos(pr_settings_t *settings, tmp_stack_t *stack)
{
	bam1_t *b = bam_init1();
	if(sam_read1(settings->in, settings->hdr, b) < 0) {
		fprintf(stderr, "Failed to read first record in bam file. Abort!\n");
		exit(EXIT_FAILURE);
	}
	while(b->core.flag & (BAM_FSUPPLEMENTARY | BAM_FSECONDARY)) {
		sam_write1(settings->out, settings->hdr, b);
		if(sam_read1(settings->in, settings->hdr, b) < 0)
		fprintf(stderr, "Failed to read first non-secondary/supplementary in bam file. Abort!\n");
		exit(EXIT_FAILURE);
	}
	stack_insert(stack, b);
	if(!(settings->in && settings->hdr)) {
		fprintf(stderr, "Failed to open input bam... WTF?\n");
		exit(EXIT_FAILURE);
	}
    while (sam_read1(settings->in, settings->hdr, b) >= 0) {
    	if(b->core.flag & (BAM_FSUPPLEMENTARY | BAM_FSECONDARY)) {
    		sam_write1(settings->out, settings->hdr, b);
    		continue;
    	}
        if(same_stack_pos(b, *stack->a)) {
        	stack_insert(stack, b);
        	continue;
        }
        else {
        	flatten_stack_linear(stack, settings); // Change this later if the chemistry necessitates it.
        	write_stack(stack, settings);
        	stack->n = 1;
        	stack->a[0] = bam_dup1(b);
        }
    }
    bam_destroy1(b);
}

static inline void pr_loop_ucs(pr_settings_t *settings, tmp_stack_t *stack)
{
#if !NDEBUG
	fprintf(stderr, "Now beginning pr_loop_ucs.\n");
#endif
	bam1_t *b = bam_init1();
	if(sam_read1(settings->in, settings->hdr, b) < 0) {
		fprintf(stderr, "Failed to read first record in bam file. Abort!\n");
		exit(EXIT_FAILURE);
	}
	while(b->core.flag & (BAM_FSUPPLEMENTARY | BAM_FSECONDARY)) {
		sam_write1(settings->out, settings->hdr, b);
		if(sam_read1(settings->in, settings->hdr, b) < 0)
		fprintf(stderr, "Failed to read first non-secondary/supplementary in bam file. Abort!\n");
		exit(EXIT_FAILURE);
	}
	stack_insert(stack, b);
	if(!(settings->in && settings->hdr)) {
		fprintf(stderr, "Failed to open input bam... WTF?\n");
		exit(EXIT_FAILURE);
	}
    while (sam_read1(settings->in, settings->hdr, b) >= 0) {
    	if(b->core.flag & (BAM_FSUPPLEMENTARY | BAM_FSECONDARY)) {
    		sam_write1(settings->out, settings->hdr, b);
    		continue;
    	}
        if(same_stack_ucs(b, *stack->a)) {
        	stack_insert(stack, b);
        	continue;
        }
        else {
        	flatten_stack_linear(stack, settings); // Change this later if the chemistry necessitates it.
        	write_stack(stack, settings);
        	stack->n = 1;
        	stack->a[0] = bam_dup1(b);
        }
    }
    bam_destroy1(b);
}


static inline void pr_loop(pr_settings_t *settings, tmp_stack_t *stack)
{
	settings->cmpkey ? pr_loop_ucs(settings, stack): pr_loop_pos(settings, stack);
}


void bam_rsq_core(pr_settings_t *settings)
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
    pr_loop(settings, &stack); // Core
    free(stack.a);
}

int pr_usage(void)
{
    fprintf(stderr, "\n");
    fprintf(stderr, "Usage:  bam_rsq [-srtu] [-f <to_realign.fq>] <input.srt.bam> <output.bam>\n\n");
    fprintf(stderr, "Option: -s    pr for SE reads [Not implemented]\n");
    fprintf(stderr, "        -r    Realign reads with no changed bases. Default: False.\n");
    fprintf(stderr, "        -t    Mismatch limit. Default: 2\n");
    fprintf(stderr, "        -l    Set bam compression level. Valid: 0-9. (0 == uncompresed)\n");
    fprintf(stderr, "        -u    Flag to use unclipped start positions instead of pos/mpos for identifying potential duplicates.\n"
    		"Note: This requires pre-processing with mark_unclipped from ppbwa (http://github.com/noseatbelts/ppbwa).\n");

    sam_global_opt_help(stderr, "-....");
    return 1;
}


int bam_rsq(int argc, char *argv[])
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
    settings->cmpkey = POS;
    settings->is_se = 0;
    settings->realign_unchanged = 0;

    char fqname[200] = "";

    while ((c = getopt_long(argc, argv, "l:f:t:aur?h", lopts, NULL)) >= 0) {
        switch (c) {
        case 'r': settings->realign_unchanged = 1; break;
        case 'u': settings->cmpkey = UCS; break;
        case 't': settings->mmlim = atoi(optarg); break;
        case 'f': strcpy(fqname, optarg); break;
        case 'l': wmode[2] = atoi(optarg) + '0'; if(wmode[2] > '9') wmode[2] = '9'; break;
        default:  if (parse_sam_global_opt(c, optarg, lopts, &ga) == 0) break;
            /* else fall-through */
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
    else bam_rsq_core(settings);
    bam_hdr_destroy(settings->hdr);
    sam_close(settings->in); sam_close(settings->out);
    if(settings->fqh)
    	fclose(settings->fqh);
    return 0;
}
