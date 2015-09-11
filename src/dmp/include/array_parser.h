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
    FILE *fp = fopen(fname, "r");
    if(!fp) {
        fprintf(stderr, "Could not open file %s. Abort mission!\n", fname);
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
    int readlen = count_lines(qual_rescale_fname);
    FILE *fp = fopen(qual_rescale_fname, "r");
    if(!fp) {
        fprintf(stderr, "Could not open file %s. Abort mission!\n", qual_rescale_fname);
        exit(EXIT_FAILURE);
    }
    char ***ret = (char ***)malloc(readlen * sizeof(char **));
    char *line = NULL;
    size_t len = 0;
    ssize_t line_length;
    int lineno = 0;
    char *qscore_tok;
    char *basecall_tok;
    size_t line_count = 0;
    size_t qs_count = 0;
    size_t bc_count = 0;
    char **qb_arr = (char **)malloc(39 * sizeof(char *));
    for(int i = 0; i < 39; ++i) {
        qb_arr[i] = (char *)malloc(4 * sizeof(char));
    }
    while ((line_length = getline(&line, &len, fp)) != -1) {
        qscore_tok = strtok(line, "\t");
        while(qscore_tok) {
            basecall_tok = strtok(qscore_tok, "\t");
            while(basecall_tok) {
                qb_arr[qs_count][bc_count] = atoi(basecall_tok);
                basecall_tok = strtok(qscore_tok, "|");
                ++bc_count;
            }
            qscore_tok = strtok(line, "\t");
            ++qs_count;
        }
        ret[line_count] = qb_arr;
        ++line_count;
    }
    if(line) free(line);
    fclose(fp);
    return ret;
}
