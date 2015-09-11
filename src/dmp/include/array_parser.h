#include <math.h>
#include <omp.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>


// Shamelessly stolen from the source code for 'wc'.
typedef unsigned long count_t; /* counter type */

inline int count_lines(char *fname) {
	int ret = 0;
	FILE *fp = fopen(fname);
	if(!fp) {
		fprintf(stdrr, "Could not open file %s. Abort mission!\n", qual_rescale_fname);
		exit(EXIT_FAILURE);
	}
	char c;
	while ((c = getc(fp)) != EOF) {
		if((c) == '\n') {
			++ret;
		}
	}
	fclose(fp);
	return ret;
}

inline char ***parse_rescaler(char *qual_rescale_fname)
{
	int n_cycles = count_lines(qual_rescale_fname);
	if(!fp) {
		fprintf(stdrr, "Could not open file %s. Abort mission!\n", qual_rescale_fname);
		exit(EXIT_FAILURE);
	}
	char ***ret = (char ***)malloc(n_cycles * sizeof(char **));
	FILE *fp = fopen(qual_rescale_fname, "r");
	char *line = NULL;
	size_t len = 0;
	ssize_t line_length;
	int lineno = 0;
	while ((line_length = getline(&line, &len, fp)) != -1) {
         //
	}
	if(line) free(line);
    fclose(fp);
	return rescaler;
}
