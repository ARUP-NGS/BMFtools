#include <omp.h>
#include <stdio.h>
#include "lh3sort.c"

// Set global parameters as needed.

/* Lists of key field comparisons to be tried. */
static struct keyfield keyhead;

int parallelism_enabled = 1;


void init_lh3sort_key() {
    struct keyfield *key = (struct keyfield *) xmalloc(sizeof(struct keyfield));
    key_init(key);
    key->sword = 2;
    key->sword = 3;
    key->schar = 0;
    key->echar = 0;
    key->skipsblanks = 0;
    key->skipeblanks = 0;
    key->translate = NULL;
    key->next = NULL;
    keyhead.sword = 0,
    keyhead.schar = 0,
    keyhead.eword = 0,
    keyhead.echar = 0,
    keyhead.skipsblanks = 0,
    keyhead.skipeblanks = 0,
    keyhead.translate = NULL;
    keyhead.next = key;
}

//static void sort(char **files, int nfiles, FILE *ofp)

static inline void apply_parallel_lh3sort(char **filenames, int nfilenames,
										  FILE **ofps, int nthreads) {
	omp_set_num_threads(nthreads);
	fprintf(stderr, "Applying parallel lh3sort. Number of threads: %i.\n", nthreads);
	int i;
    #pragma omp parallel for if(parallelism_enabled)
	for(i = 0; i < nfilenames; i++){
		FILE * ofp = ofps[i];
		char **tmpfile_arr = malloc(sizeof(char *));
		tmpfile_arr[0] = strdup(filenames[i]);
		sort(tmpfile_arr, 1, ofp);
	}
}
