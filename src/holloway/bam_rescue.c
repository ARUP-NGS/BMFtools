/*  bam_rescue.c -- duplicate read detection.

    Copyright (C) 2009, 2015 Genome Research Ltd.
    Portions copyright (C) 2009 Broad Institute.

    Author: Heng Li <lh3@sanger.ac.uk>

    Modified 9/15 by Daniel baker <daniel.baker@aruplab.com>, <d.nephi.baker@gmail.com>

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

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <zlib.h>
#include <unistd.h>
#include "htslib/sam.h"
#include "sam_opts.h"
#include "bam.h" // for bam_get_library
#include "kingfisher.h" // For igamc/igamc_pvalues
#include "cstr_utils.h"
#include "sam_opts.h"

typedef bam1_t *bam1_p;


typedef struct tmp_stack {
    int n, max;
    bam1_t **a;
} tmp_stack_t;

typedef struct rescue_settings {
    FILE *fp;
    int write_nc2bam; // Write the reads with the number of nucleotides changed to bam rather than to a fastq file.
    int hd_thresh; // Threshold for hamming distance. If < hd_thresh, reads are considered to be from the same family.
    int is_se; // Is single end
    int force_se; // Is single end
    char *fq_fname;
} rescue_settings_t;


// arr1 and arr2 should each be a uint64_t buf[2] or a malloc'd char * of size 2.
static inline int arr_cmp(char *arr1, char *arr2) {
    return (arr1[0] == arr2[0] && arr1[1] == arr2[1]);
}

#ifndef IS_REVERSE
#define IS_REVERSE(bam) ((bam)->core.flag&BAM_FREVERSE)
#endif

#ifndef IS_MATE_REVERSE
#define IS_MATE_REVERSE(bam) ((bam)->core.flag&BAM_FMREVERSE)
#endif

#ifndef IS_READ2
#define IS_READ2(bam) ((bam)->core.flag&BAM_FREAD2)
#endif

#ifndef IS_READ1
#define IS_READ1(bam) ((bam)->core.flag&BAM_FREAD1)
#endif

#ifndef ARR_CMP
#define ARR_CMP(arr1, arr2) (arr1[0] == arr2[0] && arr1[1] == arr2[1])
#endif

// FUNCTION DEFINITIONS
char *trim_ext(char *fname);
uint8_t *bam_aux_get(const bam1_t *b, const char tag[2]);
int pvalue_to_phred(double pvalue);


/*
// buf should be a uint64_t buf[2] or a malloc'd char * of size 2.
#ifndef ARR_SETKEY
#define ARR_SETKEY(bam, buf) buf = BAM2KEY(bam);
#endif
*/

static inline void ARR_SETKEY(bam1_t *bam, uint64_t buf[2])
{
    buf[0] = (((uint64_t)IS_REVERSE(bam)) << 59 | ((uint64_t)IS_MATE_REVERSE(bam)) << 57|
              ((uint64_t)IS_READ1(bam)) << 55 | ((uint64_t)IS_READ2(bam)) << 53 |
              ((uint64_t)bam->core.mtid) << 44 | ((uint64_t)bam->core.tid) << 28);
    buf[1] = (((uint64_t)(bam->core.pos) << 32 | bam->core.mpos));
    return;
}


static inline void stack_insert(tmp_stack_t *stack, bam1_t *b)
{
    if (stack->n == stack->max) {
        stack->max = stack->max? stack->max<<1 : 0x10000;
        stack->a = (bam1_t**)realloc(stack->a, sizeof(bam1_t*) * stack->max);
    }
    stack->a[stack->n++] = b;
}

    /*
     * TODO:
     * Expand the rescaling to 4 dimensions (Done.)
     * Add in commandline option key for users (Done.)
     *
     * Yet TODO
     * Fill in write_bam1_nc (for reads to be realigned)
     * Add in the count for RC to crms
     * Expand the max_phreds in all of the marksplits to make sure it's the highest quality score for the base that was called as consensus.
     */

static inline void write_bam1_nc(bam1_t *b, FILE *fp)
{
	/*
	 * TODO HERE:
	 * PV
	 * FA
	 * FM
	 * NN
	 * NC
	 * RC
	 * FP
	 */
    fprintf(stderr, "Not implemented (write_bam1_nc). Abort mission!\n");
    exit(EXIT_FAILURE);
}


static inline void write_stack(tmp_stack_t *stack, samFile *out, bam_hdr_t *hdr, rescue_settings_t *settings)
{
    int i;
    for (i = 0; i != stack->n; ++i) {
        if(stack->a[i]) { // Since I'll be clearing out the values after updating, I am checking to see if it's not NULL first.
            if(bam_aux2i(bam_aux_get(stack->a[i], "NC")) || bam_aux2i(bam_aux_get(stack->a[i], "NN"))) { // If NN or NC are non-zero.
                if(settings->write_nc2bam) {
                    sam_write1(out, hdr, stack->a[i]);
                }
                else {
                    write_bam1_nc(stack->a[i], settings->fp);
                }
            }
            else {
                sam_write1(out, hdr, stack->a[i]);
            }
            bam_destroy1(stack->a[i]);
            stack->a[i] = NULL;
        }
    }
    stack->n = 0;
}

static inline int hamming_dist_test(char *bs1, char *bs2, int hd_thresh)
{
    int mm = 0;
    for(int i = 0; bs1[i]; ++i) { // Gives up once it reaches a null terminus. Convenient, since string tags are null-terminated in the sam file format.
        if(bs1[i] != bs2[i]) {
            ++mm;
            if(!(mm - hd_thresh)) {
                return 0;
            }
        }
    }
    return 1;
}

static inline int *bam_aux2ip(uint8_t *ptr)
{
    return (int *)(++ptr);
}

static inline int disc_pvalues(double pv_better, double pv_worse)
{
    return pvalue_to_phred(igamc(2., LOG10_TO_CHI2(pv_better - (10. * log10(1 - pow(10., (pv_worse * 0.1)))))));
}

static inline void update_int_tag(bam1_t *b, const char key[2], int increment)
{
    uint8_t *tag_ptr = bam_aux_get(b, key);
    if(!tag_ptr) {
        bam_aux_append(b, key, 'i', sizeof(int), (uint8_t *)&increment);
    }
    else {
        *(int *)++tag_ptr += increment;
    }
    return;
}

static inline void update_int_ptr(uint8_t *ptr1, uint8_t *inc_ptr)
{
#if !NDEBUG
    if(!inc_ptr) {
        fprintf(stderr, "Missing RC tag. Abort mission!\n");
        exit(EXIT_FAILURE);
    }
    else {
#endif
        *(int *)(++ptr1) += *(int *)(++inc_ptr);
#if !NDEBUG
    }
#endif
    return;
}

static inline void update_bam1(bam1_t *p, bam1_t *b, FILE *fp)
{
    int n_changed = 0;
    int mask = 0;
    int *bPV = (int *)bam_aux_get(b, "PV"); // Length of this should be b->l_qseq
    int *pPV = (int *)bam_aux_get(p, "PV"); // Length of this should be b->l_qseq
    int *bFA = (int *)bam_aux_get(b, "FA"); // Length of this should be b->l_qseq
    int *pFA = (int *)bam_aux_get(p, "FA"); // Length of this should be b->l_qseq
    update_int_ptr(bam_aux_get(p, "FM"), bam_aux_get(b, "FM")); // p.FM += b.FM
#if !NDEBUG
    if(!bPV || !pPV) {
        fprintf(stderr, "Required PV tag not found. Abort mission!\n");
        exit(EXIT_FAILURE);
    }
    if(!bFA || !pFA) {
        fprintf(stderr, "Required FA tag not found. Abort mission!\n");
        exit(EXIT_FAILURE);
    }
#endif

    char *bSeq = (char *)bam_get_seq(b);
    char *pSeq = (char *)bam_get_seq(p);
    char *bQual = (char *)bam_get_qual(b);
    char *pQual = (char *)bam_get_qual(p);
    for(int i = 0; i < p->core.l_qseq; ++i) {
        if(bSeq[i] == pSeq[i]) {
            pPV[i] = igamc(2., LOG10_TO_CHI2(pPV[i] + bPV[i]) / 2.0);
            pFA[i] += bFA[i];
            if(bQual[i] > pQual[i])
                pQual[i] = bQual[i];
        }
        else {
            if(pPV[i] > bPV[i]) {
                pPV[i] = disc_pvalues(pPV[i], bPV[i]);
                // Leave pFA alone since it didn't change.
            }
            else {
                pPV[i] = disc_pvalues(bPV[i], pPV[i]);
                pSeq[i] = bSeq[i];
                pFA[i] = bFA[i];
                pQual[i] = bQual[i];
                ++n_changed;
            }
            if(pPV[i] < 3) {
                pSeq[i] = 'N';
                pFA[i] = 0;
                pPV[i] = 0;
                pQual[i] = 33; // 33 phred == 0 quality == '!'
                ++mask;
            }
        }

        // Add if not present.
        update_int_tag(p, "NC", n_changed);
        update_int_tag(p, "NN", mask);
        uint8_t *p_rc_ptr = bam_aux_get(p, "RC");
        if(p_rc_ptr) {
            update_int_ptr(p_rc_ptr, bam_aux_get(b, "RC"));
        }
        bam_destroy1(b);
        b = NULL;
    }
}

static inline void flatten_stack(tmp_stack_t *stack, rescue_settings_t *settings_ptr)
{
    uint8_t *cp, *cb;
    int pass;
    bam1_t *b, *p;
    for(int i = 0; i < stack->n; ++i) {
        b = stack->a[i];
        cb = bam_aux_get(b, "BS");
        for(int j = i + 1; j < stack->n; ++i) {
            p = stack->a[j];
            cp = bam_aux_get(p, "BS");
        }
    }
    if(hamming_dist_test(++cp, ++cb, settings_ptr->hd_thresh)) { // Increment these pointers to get to t
        update_bam1(p, b, settings_ptr->fp); // Update record p with b
    }
    return;
}

void bam_rescue_core(samFile *in, bam_hdr_t *hdr, samFile *out, rescue_settings_t *settings_ptr)
{
    bam1_t *b, *tmpb;
    int last_tid = -1;
    tmp_stack_t stack;
    uint64_t last_arr[2] = {0, 0};
    uint64_t current_arr[2] = {0, 0};

    b = bam_init1();
    memset(&stack, 0, sizeof(tmp_stack_t));

    while (sam_read1(in, hdr, b) >= 0) {
        bam1_core_t *c = &b->core;
        ARR_SETKEY(b, current_arr);
        if (!ARR_CMP(last_arr, current_arr)) {
            flatten_stack(&stack, settings_ptr);
            write_stack(&stack, out, hdr, settings_ptr); // write the result
            if ((int)c->tid < 0) { // append unmapped reads
                sam_write1(out, hdr, b);
                while (sam_read1(in, hdr, b) >= 0) sam_write1(out, hdr, b);
                break;
            }
            if(last_tid != c->tid) {
                last_tid = c->tid;
                fprintf(stderr, "[bam_rescue_core] processing reference %s...\n", hdr->target_name[c->tid]);
            }
        }
        if (!(c->flag&BAM_FPAIRED)) {
            sam_write1(out, hdr, b);
        } else { // paired, head
            tmpb = bam_dup1(b);
            stack_insert(&stack, tmpb);
        }
    }
    free(stack.a);
    bam_destroy1(b);
}

void bam_rescue_se_core(samFile *in, bam_hdr_t *hdr, samFile *out, rescue_settings_t *settings_ptr)
{
    fprintf(stderr, "Single-end rescue not implemented. Abort mission! (Or you could code it -- let me know!)\n");
    exit(EXIT_FAILURE);
}

static int rescue_usage(void) {
    fprintf(stderr, "\n");
    fprintf(stderr, "Usage:  samtools rescue [-sS] <input.srt.bam> <output.bam>\n\n");
    fprintf(stderr, "Option: -s    rescue for SE reads\n");
    fprintf(stderr, "        -S    treat PE reads as SE in rescue (force -s)\n");
    fprintf(stderr, "        -p    flag - if set, write updated bam records without realigning.\n");
    fprintf(stderr, "        -t    threshold - number of mismatches to consider barcodes distinct. [int] Default: 2\n");
    fprintf(stderr, "        -f    Fastq output path. Set to '-' for stdout. If unset, -p flag must be set to true. [int]\n");

    sam_global_opt_help(stderr, "-....");
    return 1;
}

int bam_rescue(int argc, char *argv[])
{
    /*
     *
typedef struct rescue_settings {
    int write_nc2bam; // Write the reads with the number of nucleotides changed to bam rather than to a fastq file.
    int hd_thresh; // Threshold for hamming distance. If < hd_thresh, reads are considered to be from the same family.
    int is_se; // Is single end
    int force_se; // Is single end
} rescue_settings_t;
     */
    rescue_settings_t settings = {
            .write_nc2bam = 0,
            .hd_thresh = 2,
            .force_se = 0,
            .is_se = 0,
            .fp = NULL,
            .fq_fname = NULL
    };
    samFile *in, *out;
    bam_hdr_t *header;
    char wmode[3] = {'w', 'b', 0};
    sam_global_args ga = SAM_GLOBAL_ARGS_INIT;

    static const struct option lopts[] = {
        SAM_OPT_GLOBAL_OPTIONS('-', 0, 0, 0, 0),
        { NULL, 0, NULL, 0 }
    };

    int c;
    while ((c = getopt_long(argc, argv, "sSpt:", lopts, NULL)) >= 0) {
        switch (c) {
        case 's': settings.is_se = 1; break;
        case 'S': settings.force_se = settings.is_se = 1; break;
        case 't': settings.hd_thresh = atoi(optarg); break;
        case 'p': settings.write_nc2bam = 1; break;
        case 'f': if(optarg[0] == '-') {
                      settings.fq_fname = strdup("stdout");
                      settings.fp = stdout;
                  }
                  else {
                      settings.fq_fname = trim_ext(optarg);
                      settings.fq_fname = realloc(settings.fq_fname, strlen(settings.fq_fname) + 10);
                      strcat(settings.fq_fname, ".ra.fastq");
                      settings.fp = fopen(settings.fq_fname, "w");
                  }
                  break;
        default:  if (parse_sam_global_opt(c, optarg, lopts, &ga) == 0) break;
            /* else fall-through */
        case '?': return rescue_usage();
        }
    }

    if(settings.write_nc2bam && settings.fp) {
        fprintf(stderr, "Conflicting options - write to bam and to fastq. "
                        "This will result in non-unique primary alignments "
                        "for a given read name. Abort mission!\n");
        exit(EXIT_FAILURE);
    }
    if (optind + 2 > argc)
        return rescue_usage();

    in = sam_open_format(argv[optind], "r", &ga.in);
    header = sam_hdr_read(in);
    if (header == NULL || header->n_targets == 0) {
        fprintf(stderr, "[bam_rescue] input SAM does not have header. Abort!\n");
        return 1;
    }

    sam_open_mode(wmode+1, argv[optind+1], NULL);
    out = sam_open_format(argv[optind+1], wmode, &ga.out);
    if (in == 0 || out == 0) {
        fprintf(stderr, "[bam_rescue] fail to read/write input files\n");
        return 1;
    }
    sam_hdr_write(out, header);

    if (settings.is_se) bam_rescue_se_core(in, header, out, &settings);
    else bam_rescue_core(in, header, out, &settings);
    bam_hdr_destroy(header);
    sam_close(in); sam_close(out);
    if(settings.fp) {
        fclose(settings.fp);
    };
    if(settings.fq_fname) {
        free(settings.fq_fname);
        settings.fq_fname = NULL;
    }
    return 0;
}


int main(int argc, char *argv[]) {
    return bam_rescue(argc, argv);
}
