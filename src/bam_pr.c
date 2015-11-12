/*  bam_pr.c -- duplicate read detection.

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
#include "bam_pr.h"
static inline void update_bam1(bam1_t *p, bam1_t *b)
{
	uint8_t *bdata, *pdata;
	int n_changed = 0;
	int mask = 0;
	int ra_val = 1;
	if(!b || !p) {
		// If the
		fprintf(stderr, "One of these records is null. Abort!\n");
		exit(EXIT_FAILURE);
	}
	int bFM = bam_aux2i(bam_aux_get(b, "FM"));
	int pFM = bam_aux2i(bam_aux_get(p, "FM"));
	pFM += bFM;
	bam_aux_del(p, bam_aux_get(p, "FM"));
	bam_aux_append(p, "FM", 'i', sizeof(int), (uint8_t *)&pFM);
	pdata = bam_aux_get(p, "FM");
	if(pFM < bFM)
		memcpy(bam_get_qname(p), bam_get_qname(b), b->core.l_qname);
	int pRV;
	pdata = bam_aux_get(p, "RV");
	if(pdata) {
		pRV = bam_aux2i(pdata);
		bdata = bam_aux_get(b, "RV");
		pRV += bam_aux2i(bdata);
		bam_aux_del(p, bam_aux_get(p, "RV"));
		bam_aux_append(p, "RV", 'i', sizeof(int), (uint8_t *)&pRV);
	}
	int bFP = bam_aux2i(bam_aux_get(b, "FP"));
	int pFP = bam_aux2i(bam_aux_get(p, "FP"));
	if(bFP && !pFP) {
		pFP = 1;
		bam_aux_del(p, bam_aux_get(p, "FP"));
		bam_aux_append(p, "FP", 'i', sizeof(int), (uint8_t *)&pFP);
	}
	uint32_t *bPV = (uint32_t *)array_tag(b, "PV"); // Length of this should be b->l_qseq
	uint32_t *pPV = (uint32_t *)array_tag(p, "PV"); // Length of this should be b->l_qseq
	uint32_t *bFA = (uint32_t *)array_tag(b, "FA"); // Length of this should be b->l_qseq
	uint32_t *pFA = (uint32_t *)array_tag(p, "FA"); // Length of this should be b->l_qseq
	if(!bPV || !pPV) {
		fprintf(stderr, "Required PV tag not found. Abort mission! Read names: %s, %s.\n", bam_get_qname(b), bam_get_qname(p));
		exit(EXIT_FAILURE);
	}
	if(!bFA || !pFA) {
		fprintf(stderr, "Required FA tag not found. Abort mission!\n");
		exit(EXIT_FAILURE);
	}

	char *bSeq = (char *)bam_get_seq(b);
	char *pSeq = (char *)bam_get_seq(p);
	char *bQual = (char *)bam_get_qual(b);
	char *pQual = (char *)bam_get_qual(p);
	if(!(bSeq && pSeq && bQual &&pQual)) {
		fprintf(stderr, "Qual strings or sequence strings are null. Abort!\n");
	}
	for(int i = 0; i < p->core.l_qseq; ++i) {
		if(bSeq[i] == pSeq[i]) {
			pPV[i] = (uint32_t)(-10 * log10(igamc(2., AVG_LOG_TO_CHI2(pPV[i] + bPV[i]))));
			pFA[i] += bFA[i];
			if(bQual[i] > pQual[i])
				pQual[i] = bQual[i];
		}
		else {
			if(pPV[i] > bPV[i])
				// Leave pFA alone since it didn't change.
				pPV[i] = disc_pvalues(pPV[i], bPV[i]);
			else {
				pPV[i] = disc_pvalues(bPV[i], pPV[i]);
				pSeq[i] = bSeq[i];
				pFA[i] = bFA[i];
				pQual[i] = bQual[i];
				++n_changed;
			}
		}
		if(pPV[i] < 3) {
			pSeq[i] = 'N';
			pFA[i] = 0;
			pPV[i] = 0;
			pQual[i] = '!'; // 33 phred == 0 quality == '!'
			++mask;
		}
	}
	pdata = bam_aux_get(p, "RA");
	if(!pdata) {
		bam_aux_append(p, "RA", 'i', sizeof(int), (uint8_t *)&ra_val);
		//fprintf(stderr, "Appended RA tag to read %s.\n", (char *)bam_get_qname(p));
	}
	else {
		int pRA = bam_aux2i(pdata);
		pRA |= ra_val;
		bam_aux_del(p, bam_aux_get(p, "RA"));
		bam_aux_append(p, "RA", 'i', sizeof(int), (uint8_t *)&pRA);
		//fprintf(stderr, "Deleted existing and appended RA tag to read %s.\n", (char *)bam_get_qname(p));
	}
	bam_destroy1(b);
	b = NULL;
#if !NDEBUG
	pdata = bam_aux_get(p, "RA");
	if(!pdata)
		exit(EXIT_FAILURE);
#endif
}


void bam_prse_core(pr_settings_t settings)
{
    bam1_t *b;
    int last_tid = -2;

    tmp_stack_t stack;
    resize_stack(&stack, STACK_START);

    b = bam_init1();
    while (sam_read1(settings.in, settings.hdr, b) >= 0) {
    	if(stack.n) {

    	}
    }
    bam_destroy1(b);
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

static inline int hd_annealed(bam1_t *a, bam1_t *b, int mmlim)
{
	int lhd = hd_linear(a, b, mmlim);
	return lhd ? lhd: annealed_check(a, b, mmlim);
}

static inline void flatten_stack_linear(tmp_stack_t *stack, pr_settings_t *settings)
{
	for(int i = 0; i < stack->n; ++i) {
		for(int j = i + 1; j < stack->n; ++j) {
			if(hd_linear(stack->a[i], stack->a[j])) {
				update_bam1(stack->a[j], stack->a[i]);
				break;
				// "break" in case there are multiple within hamming distance.
				// Otherwise, I'll end up having memory mistakes.
				// Besides, that read set will get merged into the later read in the set.
			}
		}
	}
}

static inline void flatten_stack_annealed(tmp_stack_t *stack, pr_settings_t *settings)
{/*
	int hd;
	for(int i = 0; i < stack->n; ++i) {
		for(int j = i + 1; j < stack->n; ++j) {
			if((hd = hd_linear(stack->a[i], stack->a[j]))) {
				if(hd > 0)
					update_bam1(stack->a[j], stack->a[i]);
				else
					update_bam1_rc(stack->a[j], stack->a[i]);
				break;
				// "break" in case there are multiple within hamming distance.
				// Otherwise, I'll end up having memory mistakes.
				// Besides, that read set will get merged into the later read in the set.
			}
		}
	}
	*/
}

static inline void flatten_stack(tmp_stack_t *stack, pr_settings_t *settings)
{
	settings->annealed ? flatten_stack_annealed(stack, settings): flatten_stack_linear(stack, settings);
}

static inline void

static inline void pr_loop_pos(pr_settings_t *settings, tmp_stack_t *stack)
{
    while (sam_read1(settings.in, settings.hdr, b) >= 0) {
        if(same_stack_pos(b, stack->a[0])) {
        	stack_insert(stack, b);
        	continue;
        }
        else {
        	flatten_stack_linear(stack, settings);
        	write_stack(stack, settings);
        	resize_stack(stack, 1);
        	stack->a[0] = bam_dup1(b);
        }
    }
}

static inline void pr_loop_ucs(pr_settings_t *settings, tmp_stack_t *stack)
{
    while (sam_read1(settings->in, settings->hdr, b) >= 0) {
        if(same_stack_ucs(b, stack->a[0])) {
        	stack_insert(stack, b);
        	continue;
        }
        else {
        	flatten_stack_linear(stack, settings); // Change this later if the chemistry necessitates it.
        	write_stack(stack, settings);
        	resize_stack(stack, 1);
        	stack->a[0] = bam_dup1(b);
        }
    }
}


static inline void pr_loop(pr_settings_t *settings, tmp_stack_t *stack)
{
	settings.cmpkey ? pr_loop_ucs(settings, stack): pr_loop_pos(settings, stack);
}


void bam_pr_core(pr_settings_t settings)
{
    bam1_t *b;
    int last_tid = -1, last_pos = -1;
    tmp_stack_t stack;

    b = bam_init1();
    resize_stack(&stack, STACK_START);
    memset(&stack, 0, sizeof(tmp_stack_t) * stack->max);
    pr_loop(&settings, &stack); // Core
    free(stack.a);
    bam_destroy1(b);
}

int pr_usage(void)
{
    fprintf(stderr, "\n");
    fprintf(stderr, "Usage:  samtools pr [-sS] <input.srt.bam> <output.bam>\n\n");
    fprintf(stderr, "Option: -s    pr for SE reads\n");
    fprintf(stderr, "        -a    Use annealed rescue\n");
    fprintf(stderr, "        -t    Mismatch limit. Default: 2\n");
    fprintf(stderr, "        -u    Flag to use unclipped start positions instead of pos/mpos for identifying potential duplicates.\n"
    		"Note: This requires pre-processing with mark_unclipped from ppbwa (http://github.com/noseatbelts/ppbwa).\n");

    sam_global_opt_help(stderr, "-....");
    return 1;
}

int main(int argc, char *argv[])
{
	return bam_pr(argc, argv);
}

int bam_pr(int argc, char *argv[])
{
    int c, is_se = 0, force_se = 0;
    char wmode[3] = {'w', 'b', 0};
    sam_global_args ga = SAM_GLOBAL_ARGS_INIT;

    static const struct option lopts[] = {
        SAM_OPT_GLOBAL_OPTIONS('-', 0, 0, 0, 0),
        { NULL, 0, NULL, 0 }
    };

    pr_settings_t settings = {
    		.fqh = NULL,
			.in = NULL,
			.out = NULL,
			.hdr = NULL,
			.mmlim = 2,
			.cmpkey = POS,
			.is_se = 0,
			.annealed = 0
    };

    while ((c = getopt_long(argc, argv, "t:au?h", lopts, NULL)) >= 0) {
        switch (c) {
        case 'a': settings.annealed = 1; break;
        case 'u': settings.cmpkey = UCS; break;
        case 't': settings.mmlim = atoi(optarg); break;
        default:  if (parse_sam_global_opt(c, optarg, lopts, &ga) == 0) break;
            /* else fall-through */
        case '?': return pr_usage();
        }
    }
    if (optind + 2 > argc)
        return pr_usage();

    settings.in = sam_open_format(argv[optind], "r", &ga.in);
    settings.hdr = sam_hdr_read(in);
    if (settings.hdr == NULL || settings.hdr->n_targets == 0) {
        fprintf(stderr, "[bam_pr] input SAM does not have header. Abort!\n");
        return 1;
    }

    sam_open_mode(wmode+1, argv[optind+1], NULL);
    settings.out = sam_open_format(argv[optind+1], wmode, &ga.out);
    if (setings.in == 0 || settings.out == 0) {
        fprintf(stderr, "[bam_pr] fail to read/write input files\n");
        return 1;
    }
    sam_hdr_write(settings.out, settings.out);

    if (is_se) bam_prse_core(settings);
    else bam_pr_core(settings);
    bam_hdr_destroy(settings.hdr);
    sam_close(settings.in); sam_close(settings.out);
    return 0;
}
