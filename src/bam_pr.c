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
#if !NDEBUG
	//fprintf(stderr, "Now updating bam1_t with name %s with bam1_t with name %s.\n", bam_get_qname(p), bam_get_qname(b));
#endif
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
	for(int i = 0; i < p->core.l_qseq; ++i) {
		if(bam_seqi(pSeq, i) == bam_seqi( bSeq, i)) {
			pPV[i] = (uint32_t)(-10 * log10(igamc(2., AVG_LOG_TO_CHI2(pPV[i], bPV[i]))));
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
				set_base(pSeq, bSeq, i);
				pFA[i] = bFA[i];
				pQual[i] = bQual[i];
				++n_changed;
			}
		}
		if(pPV[i] < 3) {
			pSeq[(i)>>1] |= (0xf << ((!(i % 2)) * 4)); // Set the base to N
			pFA[i] = 0;
			pPV[i] = 0;
			pQual[i] = 0; // Note: this is not shifted by 33.
			++mask;
		}
	}
	pdata = bam_aux_get(p, "RA");
	if(!pdata) {
		bam_aux_append(p, "RA", 'i', sizeof(int), (uint8_t *)&ra_val);
		//fprintf(stderr, "Appended RA tag to read %s.\n", (char *)bam_get_qname(p));
	}
	bam_destroy1(b);
	b = NULL;
#if !NDEBUG
	pdata = bam_aux_get(p, "RA");
	if(!pdata)
		exit(EXIT_FAILURE);
#endif
}


void bam_prse_core(pr_settings_t *settings)
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

static inline int hd_annealed(bam1_t *a, bam1_t *b, int mmlim)
{
	int lhd = hd_linear(a, b, mmlim);
	return lhd ? lhd: annealed_check(a, b, mmlim);
}

static inline void flatten_stack_linear(tmp_stack_t *stack, pr_settings_t *settings)
{
	for(int i = 0; i < stack->n; ++i) {
		for(int j = i + 1; j < stack->n; ++j) {
			if(hd_linear(stack->a[i], stack->a[j], settings->mmlim)) {
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

static inline void pr_loop_pos(pr_settings_t *settings, tmp_stack_t *stack)
{
	bam1_t *b = bam_init1();
	if(sam_read1(settings->in, settings->hdr, b) < 0) {
		fprintf(stderr, "Failed to read first record in bam file. Abort!\n");
		exit(EXIT_FAILURE);
	}
	stack_insert(stack, b);
    while (sam_read1(settings->in, settings->hdr, b) >= 0) {
    	fprintf(stderr, "Current sequence: %s\n", bam_get_seq(b));
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
    bam_destroy1(b);
}

static inline void pr_loop_ucs(pr_settings_t *settings, tmp_stack_t *stack)
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
	if(!stack->a) {
		fprintf(stderr, "Looks like my stack wasn't properly initialized.\n");
		exit(EXIT_FAILURE);
	}
    while (sam_read1(settings->in, settings->hdr, b) >= 0) {
    	if(b->core.flag & (BAM_FSUPPLEMENTARY | BAM_FSECONDARY)) {
    		sam_write1(settings->out, settings->hdr, b);
    		continue;
    	}
        if(same_stack_ucs(b, stack->a[stack->n - 1])) {
#if !NDEBUG
        	if(strcmp(bam_get_qname(b), bam_get_qname(stack->a[0])) == 0) {
        		fprintf(stderr, "We're comparing records at %p and %p which have the same name. Abort!\n", b, stack->a[0]);
        		exit(EXIT_FAILURE);
        	}
#endif
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
	pr_loop_ucs(settings, stack);
	//settings->cmpkey ? pr_loop_ucs(settings, stack): pr_loop_pos(settings, stack);
}


void bam_pr_core(pr_settings_t *settings)
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
    int c;
    char wmode[3] = {'w', 'b', 0};
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
    settings->cmpkey = 0;
    settings->is_se = 0;
    settings->annealed = 0;

    char fqname[200] = "";

    while ((c = getopt_long(argc, argv, "f:t:au?h", lopts, NULL)) >= 0) {
        switch (c) {
        case 'a': settings->annealed = 1; break;
        case 'u': settings->cmpkey = 1; break;
        case 't': settings->mmlim = atoi(optarg); break;
        case 'f': strcpy(fqname, optarg); break;
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
        fprintf(stderr, "[bam_pr] input SAM does not have header. Abort!\n");
        return 1;
    }

    sam_open_mode(wmode+1, argv[optind+1], NULL);
    settings->out = sam_open_format(argv[optind+1], wmode, &ga.out);
    if (settings->in == 0 || settings->out == 0) {
        fprintf(stderr, "[bam_pr] fail to read/write input files\n");
        return 1;
    }
    sam_hdr_write(settings->out, settings->hdr);

    if (settings->is_se) bam_prse_core(settings);
    else bam_pr_core(settings);
    bam_hdr_destroy(settings->hdr);
    sam_close(settings->in); sam_close(settings->out);
    if(settings->fqh)
    	fclose(settings->fqh);
    return 0;
}
