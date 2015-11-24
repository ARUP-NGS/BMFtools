#include "bed_util.h"


khash_t(bed) *parse_bed(char *path, bam_hdr_t *header, int padding)
{
	khash_t(bed) *ret = kh_init(bed);
	FILE *ifp = fopen(path, "r");
	char *line = NULL;
	char *tok = NULL;
	size_t len = 0;
	ssize_t read;
	uint32_t tid, start, stop;
	int khr;
	khint_t k;
	while ((read = getline(&line, &len, ifp)) != -1) {
		if(line[0] == '\0' || line[0] == '#') // Empty line or comment line
			continue;
		tok = strtok(line, "\t");
		tid = (uint32_t)bam_name2id(header, tok);
#if !NDEBUG
		fprintf(stderr, "Transcript id for tok %s is %"PRIu32"\n", tok, tid);
#endif
		tok = strtok(NULL, "\t");
		start = strtoul(tok, NULL, 10);
		tok = strtok(NULL, "\t");
		stop = strtoul(tok, NULL, 10);
		k = kh_get(bed, ret, tid);
		if(k == kh_end(ret)) {
#if !NDEBUG
			fprintf(stderr, "New contig in bed hashmap: %"PRIu32".\n", tid);
#endif
			k = kh_put(bed, ret, tid, &khr);
			kh_val(ret, k).intervals = (interval_t *)calloc(1, sizeof(interval_t));
			kh_val(ret, k).intervals[0].start = start - padding > 0 ? start - padding : 0;
			kh_val(ret, k).intervals[0].end = stop + padding;
			kh_val(ret, k).n = 1;
		}
		else {
			kh_val(ret, k).intervals = (interval_t *)realloc(kh_val(ret, k).intervals, ++kh_val(ret, k).n * sizeof(interval_t));
			if(!kh_val(ret, k).intervals) {
				fprintf(stderr, "Could not allocate memory. Abort mission!\n");
				exit(EXIT_FAILURE);
			}
			kh_val(ret, k).intervals[kh_val(ret, k).n - 1].start = start;
			kh_val(ret, k).intervals[kh_val(ret, k).n - 1].end = stop;
#if !NDEBUG
			fprintf(stderr, "Number of intervals in bed file for contig %"PRIu32": %"PRIu64"\n", tid, kh_val(ret, k).n);
#endif
		}
	}
	sort_bed(ret);
	return ret;
}


/*
 * Compare two interval objects
 */
int intcmp(const void *a, const void *b)
{
	size_t diff = (*(interval_t *)a).start - (*(interval_t *)b).start;
	return (diff) ? diff: (*(interval_t *)a).end - (*(interval_t *)b).end;
}


void bed_destroy(khash_t(bed) *b)
{
	khint_t ki;
	for(ki = kh_begin(b); ki != kh_end(b); ++ki) {
		if(!kh_exist(b, ki))
			continue;
		cond_free(kh_val(b, ki).intervals);
		kh_val(b, ki).n = 0;
	}
	cond_free(b);
}

void sort_bed(khash_t(bed) *bed)
{
	khint_t k;
	for(k = kh_begin(bed); k != kh_end(bed); ++k) {
		if(!kh_exist(bed, k))
			continue;
		qsort(kh_val(bed, k).intervals, kh_val(bed, k).n, sizeof(region_set_t), &intcmp);
	}
	return;
}

