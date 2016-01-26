/*  bmf_depth.c -- bedcov subcommand.

	Modified from bedcov.c

ORIGINAL COPYRIGHT:

	Copyright (C) 2012 Broad Institute.
	Copyright (C) 2013-2014 Genome Research Ltd.

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
DEALINGS IN THE SOFTWARE. */

#include "bmf_depth.h"

#define NO_ID_STR "L'Innommable"


typedef struct {
	htsFile *fp;
	bam_hdr_t *header;
	hts_itr_t *iter;
	uint64_t *raw_counts;
	uint64_t *dmp_counts;
	int minMQ;
	int minFM;
} aux_t;

void depth_usage(int retcode)
{
	fprintf(stderr, "Usage: bmftools depth [options] -b <in.bed> <in1.bam> [...]\n\n");
	fprintf(stderr, "  -Q INT		Only count bases of at least INT quality [0]\n");
	fprintf(stderr, "  -f INT		Only count bases of at least INT Famly size (unmarked reads have FM 1) [0]\n");
	fprintf(stderr, "  -m INT		Max depth. Default: %i.\n", DEFAULT_MAX_DEPTH);
	sam_global_opt_help(stderr, "-.--.");
	exit(retcode);
}

int plp_fm_sum(const bam_pileup1_t *stack, int n_plp)
{
	int ret = 0;
	for(; n_plp; --n_plp) ret += get_fm((*stack++).b);
	return ret;
}

static int read_bam(void *data, bam1_t *b)
{
	aux_t *aux = (aux_t*)data; // data in fact is a pointer to an auxiliary structure
	int ret;
	for(;;)
	{
		ret = aux->iter? sam_itr_next(aux->fp, aux->iter, b) : sam_read1(aux->fp, aux->header, b);
		if ( ret<0 ) break;
		uint8_t *data = bam_aux_get(b, "FM");
		if ((b->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP)) ||
			(int)b->core.qual < aux->minMQ || (data && bam_aux2i(data) < aux->minFM))
				continue;
		break;
	}
	return(ret);
}

int depth_main(int argc, char *argv[])
{
	gzFile fp;
	kstring_t str;
	kstream_t *ks;
	hts_idx_t **idx;
	aux_t **aux;
	char **col_names;
	int *n_plp, dret, i, n, c, minMQ = 0;
	uint64_t *counts;
	const bam_pileup1_t **plp;
	int usage = 0, max_depth = DEFAULT_MAX_DEPTH, minFM = 0;
	char *bedpath = NULL;

	sam_global_args ga = SAM_GLOBAL_ARGS_INIT;
	static const struct option lopts[] = {
		SAM_OPT_GLOBAL_OPTIONS('-', 0, '-', '-', 0),
		{ NULL, 0, NULL, 0 }
	};

	while ((c = getopt_long(argc, argv, "Q:b:m:f:?h", lopts, NULL)) >= 0) {
		switch (c) {
		case 'Q': minMQ = atoi(optarg); break;
		case 'b': bedpath = strdup(optarg); break;
		case 'm': max_depth = atoi(optarg); break;
		case 'f': minFM = atoi(optarg); break;
		default:  if (parse_sam_global_opt(c, optarg, lopts, &ga) == 0) break;
				  /* else fall-through */
		case 'h': /* fall-through */
		case '?': usage = 1; break;
		}
		if (usage) break;
	}
	if (usage || optind > argc) // Require at least one bam
		depth_usage(EXIT_FAILURE);
	memset(&str, 0, sizeof(kstring_t));
	n = argc - optind - 1;
	aux = calloc(n, sizeof(aux_t*));
	idx = calloc(n, sizeof(hts_idx_t*));
	for (i = 0; i < n; ++i) {
		aux[i] = calloc(1, sizeof(aux_t));
		aux[i]->minMQ = minMQ;
		aux[i]->minFM = minFM;
		aux[i]->fp = sam_open_format(argv[i + optind], "r", &ga.in);
		if (aux[i]->fp)
			idx[i] = sam_index_load(aux[i]->fp, argv[i + optind]);
		if (aux[i]->fp == 0 || idx[i] == 0) {
			fprintf(stderr, "ERROR: fail to open index BAM file '%s'\n", argv[i + optind]);
			return 2;
		}
		// TODO bgzf_set_cache_size(aux[i]->fp, 20);
		aux[i]->header = sam_hdr_read(aux[i]->fp);
		if (aux[i]->header == NULL) {
			fprintf(stderr, "ERROR: failed to read header for '%s'\n",
					argv[i+optind+1]);
			return 2;
		}
	}
	if(!bedpath) {
		LOG_ERROR("Bed path required. Abort!\n");
	}
	counts = calloc(n, sizeof(uint64_t));
	int n_cols = count_lines(bedpath);
	col_names = calloc(n_cols, sizeof(char *));

	fp = gzopen(argv[optind], "rb");
	ks = ks_init(fp);
	n_plp = calloc(n, sizeof(int));
	plp = calloc(n, sizeof(bam_pileup1_t*));
	int line_num = 0;
	while (ks_getuntil(ks, KS_SEP_LINE, &str, &dret) >= 0) {
		char *p, *q;
		int tid, beg, end, pos, region_len;
		bam_mplp_t mplp;

		for (p = q = str.s; *p && *p != '\t'; ++p);
		if (*p != '\t') goto bed_error;
		*p = 0; tid = bam_name2id(aux[0]->header, q); *p = '\t';
		if (tid < 0) goto bed_error;
		for (q = p = p + 1; isdigit(*p); ++p);
		if (*p != '\t') goto bed_error;
		*p = 0; beg = atoi(q); *p = '\t';
		for (q = p = p + 1; isdigit(*p); ++p);
		if (*p == '\t' || *p == 0) {
			int c = *p;
			*p = 0; end = atoi(q); *p = c;
		} else goto bed_error;
		region_len = end - beg;
		for(i = 0; i < n; ++i)
			crealloc(aux[i]->dmp_counts, sizeof(uint64_t) * region_len),
			crealloc(aux[i]->raw_counts, sizeof(uint64_t) * region_len);
		if(*p == '\t') {
			q = ++p;
			while(*q != '\t' && *q != '\n') ++q;
			int c = *q; *q = '\0';
			restrdup(col_names[i], p);
			*q = c;
		} else restrdup(col_names[i], NO_ID_STR);

		for (i = 0; i < n; ++i) {
			if (aux[i]->iter) hts_itr_destroy(aux[i]->iter);
			aux[i]->iter = sam_itr_queryi(idx[i], tid, beg, end);
		}
		mplp = bam_mplp_init(n, read_bam, (void**)aux);
		bam_mplp_set_maxcnt(mplp, max_depth);
		memset(counts, 0, sizeof(uint64_t) * n);
		int arr_ind = 0;
		while (bam_mplp_auto(mplp, &tid, &pos, n_plp, plp) > 0) {
			if (pos >= beg && pos < end) {
				for (i = 0; i < n; ++i) {
					counts += n_plp[i];
					aux[i]->dmp_counts[arr_ind] = n_plp[i];
					aux[i]->raw_counts[arr_ind] = plp_fm_sum(mplp->plp[i], n_plp[i]);
				}
				++arr_ind; // Increment for positions in range.
			}
		}


		kputc('\t', &str);
		kputs(col_names[i], &str);
		for (i = 0; i < n; ++i) {
			kputc('\t', &str);
			kputl(counts[i], &str);
		}
		puts(str.s);
		bam_mplp_destroy(mplp);
		++line_num;
		continue;

bed_error:
		fprintf(stderr, "Errors in BED line '%s'\n", str.s);
	}
	free(n_plp); free(plp);
	ks_destroy(ks);
	gzclose(fp);

	for (i = 0; i < n; ++i) {
		cond_free(aux[i]);
		if (aux[i]->iter) hts_itr_destroy(aux[i]->iter);
		hts_idx_destroy(idx[i]);
		bam_hdr_destroy(aux[i]->header);
		sam_close(aux[i]->fp);
	}
	for(i = 0; i < n_cols; ++i) free(col_names[i]);
	free(counts);
	free(col_names);
	free(aux); free(idx);
	free(str.s);
	sam_global_args_free(&ga);
	return(EXIT_SUCCESS);
}
