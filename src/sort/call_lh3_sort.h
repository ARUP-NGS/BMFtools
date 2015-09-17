#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include "lh3sort.c"

static inline void apply_parallel_lh3sort(char **filenames, int nfilenames,
                                          FILE **ofps, int nthreads) {
    omp_set_num_threads(nthreads);
    fprintf(stderr, "Applying parallel lh3sort. Number of threads: %i.\n", nthreads);
    int i;
    #pragma omp parallel for
    for(i = 0; i < nfilenames; i++){
        FILE * ofp = ofps[i];
        char **tmpfile_arr = malloc(sizeof(char *));
        tmpfile_arr[0] = strdup(filenames[i]);
        sort(tmpfile_arr, 1, ofp);
        free(tmpfile_arr[0]);
        free(tmpfile_arr);
        xfclose(ofp);
    }
}
