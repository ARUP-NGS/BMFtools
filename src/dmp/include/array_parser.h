#pragma once

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

inline char ****parse_rescaler(char *qual_rescale_fname)
{
	int readlen = count_lines(qual_rescale_fname);
	FILE *fp = fopen(qual_rescale_fname, "r");
	if(!fp) {
		fprintf(stderr, "Could not open file %s. Abort mission!\n", qual_rescale_fname);
		exit(EXIT_FAILURE);
	}
	char ****ret = (char ****)malloc(2 * sizeof(char ***));
	for(int i = 0; i < 2; ++i) {
		ret[i] = (char ***)malloc(readlen * sizeof(char **));
		for(int j = 0; j < readlen; ++j) {
			ret[i][j] = (char **)malloc(39 * sizeof(char *));
			for(int k = 0; k < 39; ++k) {
				ret[i][j][k] = (char *)malloc(4 * sizeof(char));
			}
		}
	}
	char *line = NULL;
	size_t len = 0;
	ssize_t line_length;
	int lineno = 0;
	char *readnum_tok;
	char *qscore_tok;
	char *basecall_tok;
	size_t readnum_count = 0;
	size_t line_count = 0;
	size_t qs_count = 0;
	size_t bc_count = 0;
	while ((line_length = getline(&line, &len, fp)) != -1) {
		readnum_tok = strtok(line, "|");
		while(readnum_tok) {
			qscore_tok = strtok(readnum_tok, ",");
			while(qscore_tok) {
				basecall_tok = strtok(qscore_tok, ":");
				while(basecall_tok) {
					ret[readnum_count][line_count][qs_count][bc_count] = atoi(basecall_tok);
					basecall_tok = strtok(NULL, ":");
					++bc_count;
				}
				qscore_tok = strtok(NULL, ",");
				++qs_count;
			}
			readnum_tok = strtok(NULL, "|");
			++readnum_count;
		}
		++line_count;
	}
	if(line) free(line);
	fclose(fp);
	return ret;
}

void period_to_null(char *instr)
{
	int i = 0;
	while(instr[i]) {
		if(instr[i] == '.') {
			instr[i] = '\0';
		}
		++i;
	}
}


inline char *parse_1d_rescaler(char *qual_rescale_fname)
{
	int readlen = count_lines(qual_rescale_fname);
#if !NDEBUG
	fprintf(stderr, "Number of lines: %i.\n", readlen);
#endif
	FILE *fp = fopen(qual_rescale_fname, "r");
	if(!fp) {
		fprintf(stderr, "Could not open file %s. Abort mission!\n", qual_rescale_fname);
		exit(EXIT_FAILURE);
	}
	int arr_len = 2 * readlen * 39 * 4 * sizeof(char);
	char *ret = (char *)malloc(2 * readlen * 39 * 4 * sizeof(char));
	char *line = NULL;
	size_t len = 0;
	ssize_t line_length = 0;
	char *tok = NULL;
	size_t index = 0;
	int lnum = 0;
	fprintf(stderr, "Parsing in array from %s...\n", qual_rescale_fname);
	while ((line_length = getline(&line, &len, fp)) != -1) {
		for(int readnum = 0; readnum < 2; ++readnum) {
			for(int qnum = 0; qnum < 39; ++qnum) {
				for(int bnum = 0; bnum < 4; ++bnum) {
					tok = strtok(tok ? NULL : line, "|:,");
					if(!tok) {
						break;
					}
					period_to_null(tok);
					ret[index++] = atoi(tok);
					//fprintf(stderr, "New qual str: %s. New qual: %c. Cycle: %i. Read num: %i. Qscore num %i. Basecall %i. .\n", tok, atoi(tok) + 33, lnum + 1, readnum + 1, qnum + 2, bnum + 1);
				}
			}
		}
		++lnum;
	}
	cond_free(line);
	fclose(fp);
	return ret;
}
