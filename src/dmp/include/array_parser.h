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

char num2nuc(int num);

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

inline char num2nuc(int num)
{
	switch(num) {
	case 0:
		return 'A';
	case 1:
		return 'C';
	case 2:
		return 'G';
	case 3:
		return 'T';
	default: return 'N';
	}
}


inline char *parse_1d_rescaler(char *qual_rescale_fname)
{
	int readlen = count_lines(qual_rescale_fname);
#if !NDEBUG
	fprintf(stderr, "Number of lines: %i.\n", readlen);
#endif
	FILE *fp = fopen(qual_rescale_fname, "rb");
	if(!fp) {
		fprintf(stderr, "Could not open file %s. Abort mission!\n", qual_rescale_fname);
		exit(EXIT_FAILURE);
	}
	fseek(fp, 0, SEEK_END);
	int length;
	length = ftell(fp);
	char *buffer;
	fseek(fp, 0, SEEK_SET);
	  buffer = malloc(length * sizeof(char));
	  if(buffer) {
	    fread(buffer, 1, length, fp);
	  }
	  else {
		  fprintf(stderr, "MEMORY ERROR.\n");
		  exit(EXIT_FAILURE);
	  }
	  fclose(fp);
	int arr_len = 2 * readlen * 39 * 4 * sizeof(char);
	char *ret = (char *)malloc(2 * readlen * 39 * 4 * sizeof(char));
	char *line = NULL;
	char *tok = NULL;
	size_t index = 0;
	fprintf(stderr, "Parsing in array from %s...\n", qual_rescale_fname);
	int lnum;
	for(lnum = 0; lnum < readlen; ++lnum) {
		fprintf(stderr, "Now reading line %i.\n", lnum);
		for(int readnum = 0; readnum < 2; ++readnum) {
			fprintf(stderr, "Now working with read %i.\n", readnum + 1);
			for(int qnum = 0; qnum < 39; ++qnum) {
				fprintf(stderr, "Now working with qscore %i.\n", qnum + 2);
				for(int bnum = 0; bnum < 4; ++bnum) {
					fprintf(stderr, "Now working with base %c.\n", num2nuc(bnum));
					tok = strtok(tok ? NULL : buffer, "|:,");
					if(!tok) {
						break;
					}
					fprintf(stderr, "1d rescaler starting. (%i, %s).\n", atoi(tok), tok);
					period_to_null(tok);
					ret[index++] = (char)((int)atof(tok));
					if(ret[index - 1] < 0) {
						fprintf(stderr, "Negative integer in 1d rescaler (%i, %s). Index: %i. Abort mission!\n", ret[index - 1], tok, index - 1);
						exit(EXIT_FAILURE);
					}
					fprintf(stderr, "%i,%i,%i,%c", lnum, readnum + 1, qnum + 2, num2nuc(bnum));
					//fprintf(stderr, "New qual str: %s. New qual: %c. Cycle: %i. Read num: %i. Qscore num %i. Basecall %i. .\n", tok, atoi(tok) + 33, lnum + 1, readnum + 1, qnum + 2, bnum + 1);
				}
			}
		}
	}
	fprintf(stderr, "Final index: %i. Expected: %i.\n", index, arr_len - 1);
	fprintf(stderr, "Final line num: %i. Expected: %i.\n", lnum, readlen - 1);
	cond_free(line);
	fclose(fp);
	return ret;
}
