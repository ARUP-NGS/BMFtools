#ifndef ARRAY_PARSER_H
#define ARRAY_PARSER_H
//#include <stdlib.h>
#include <string.h>
#include "dlib/char_util.h"
#include "dlib/io_util.h"

#define nqscores 39uL // Number of q scores in sequencing.


#define free_rescaler_array(settings) \
        int _readlen = count_lines(settings.rescaler_path);\
        for(int i##settings = 0; i##settings < 2; ++i##settings) {\
            for(int j##settings = 0; j##settings < _readlen; ++j##settings) {\
                for(unsigned k##settings = 0; k##settings < nqscores; ++k##settings) {\
                    cond_free(settings.rescaler[i##settings][j##settings][k##settings]);\
                }\
                cond_free(settings.rescaler[i##settings][j##settings]);\
            }\
            cond_free(settings.rescaler[i##settings]);\
        }\
        cond_free(settings.rescaler)\

/*
 * @func rescale_qscore
 * Returns a rescaled qscore using readnumber, qscore, cycle, base call, and an associated rescaling array.
 * :param: readnum [int] 0 if 1, 1 if read 2
 * :param: qscore [char] Character in quality string, with 33 offset. (e.g., '#' means 0.)
 * :param: cycle [int] Cycle number on the sequencer, 0-based.
 * :param: base [char] base call at position.
 * @discussion: Note - the rescaled q score in the flat text file has not been offset by 33.
 * This was done to make it more human-readable, though it does necessitate incrementing the quality score
 * by 33 each time a qscore is returned. Thoughts?
 *
 */
CONST static inline char rescale_qscore(int readnum, char qscore, int cycle, char base, int readlen, char *rescaler)
{
    int index = readnum;
    int mult = 2;
    //fprintf(stderr, "index value is now: %i, mult %i.\n", index, mult);
    index += cycle * mult;
    mult *= readlen;
    //fprintf(stderr, "index value is now: %i, mult %i. Qscore: %i, Qscore index%i.\n", index, mult, qscore, qscore - 35);
    index += (qscore - 35) * mult; // Subtract 35 - 33 to get to phred space, 2 to offset by 2.
    mult *= nqscores;
    //fprintf(stderr, "index value is now: %i, mult %i.\n", index, mult);
    index += mult * nuc2num_acgt(base);
    //fprintf(stderr, "Index = %i.\n", index);
#if DBG
    if(index >= readlen * 2 * nqscores * 4 || index < 0) {
        fprintf(stderr, "Something's wrong. Index (%i) is too big or negative! Max: %i.\n", index, readlen * 2 * nqscores * 4);
        exit(EXIT_FAILURE);
    }
    //fprintf(stderr, "Value at index: %i (%c).\n", rescaler[index], rescaler[index] + 33);
    if(rescaler[index] < 0) {
        fprintf(stderr, "WTF THIS CAN'T BE BELOW 0 (%i).\n", rescaler[index]);
        exit(EXIT_FAILURE);
    }
#endif
    return rescaler[index] + 33;
}


// Shamelessly stolen from the source code for 'wc'.
typedef unsigned long count_t; /* counter type */

static inline char ****parse_rescaler(char *qual_rescale_fname)
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
            ret[i][j] = (char **)malloc(nqscores * sizeof(char *));
            for(unsigned k = 0; k < nqscores; ++k) {
                ret[i][j][k] = (char *)malloc(4 * sizeof(char));
            }
        }
    }
    char *line = NULL;
    size_t len = 0;
    ssize_t line_length;
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

static void period_to_null(char *instr)
{
    /*
    start:
    switch(*instr++) {
    case '\0':
        return;
    case '.':
        *(instr - 1) = '\0';
        return;
    }
    goto start;
    */
    while(*instr) {
        if(*(instr) == '.')
            *(instr) = '\0';
        ++instr;
    }
}


static inline char *parse_1d_rescaler(char *qual_rescale_fname)
{
    const int readlen = count_lines(qual_rescale_fname);
#if DBG
    fprintf(stderr, "[D:%s] Number of lines: %i.\n", __func__, readlen);
#endif
    FILE *fp = fopen(qual_rescale_fname, "rb");
    if(!fp) {
        fprintf(stderr, "Could not open file %s. Abort mission!\n", qual_rescale_fname);
        exit(EXIT_FAILURE);
    }
    fseek(fp, 0, SEEK_END);
    const int length = ftell(fp);
    fseek(fp, 0, SEEK_SET);
    char *buffer = (char *)malloc(length * sizeof(char));
    if(buffer)
          fread(buffer, 1, length, fp);
    else {
          fprintf(stderr, "MEMORY ERROR.\n");
          exit(EXIT_FAILURE);
    }
    fclose(fp), fp = NULL;
    const size_t arr_len = 2 * readlen * nqscores * 4;
    char *ret = (char *)malloc(arr_len * sizeof(char));
    memset(ret, -127, arr_len); // Set all of these char values to -127, which is definitely unprintable
    char *tok = NULL;
    size_t index = 0;
    fprintf(stderr, "[D:%s] Parsing in array with read len %lu from %s...\n", __func__,
            arr_len, qual_rescale_fname);
    int lnum;
    for(lnum = 0; lnum < readlen; ++lnum) {
        for(int readnum = 0; readnum < 2; ++readnum) {
            for(unsigned qnum = 0; qnum < nqscores; ++qnum) {
                for(int bnum = 0; bnum < 4; ++bnum) {
                    tok = strtok(tok ? NULL : buffer, "|:,\n");
                    if(!tok)
                        break;
                    period_to_null(tok);
                    ret[index++] = (char)((int)(atof(tok) + 0.5)); // Round up.
                    if(ret[index - 1] < 0) {
                        fprintf(stderr,
                                "[E:%s] Negative integer in 1d rescaler %i, %s). Index:"
                                " %lu. Abort mission!\n", __func__, ret[index - 1], tok, index - 1);
                        exit(EXIT_FAILURE);
                    } else if(ret[index - 1] > 93) {
                        fprintf(stderr,
                                "[W:%s] rescaled quality score above the max"
                                " that can be held in the fastq format. Capping at 93. Value: %i.\n",
                                __func__, ret[index - 1]);
                        ret[index - 1] = 93;
                    } else if(ret[index - 1] < 2) ret[index - 1] = 2;
                }
            }
        }
    }
    if(index != arr_len) {
        fprintf(stderr, "[E:%s] Unexpected index! Abort!\n", __func__);
        exit(EXIT_FAILURE);
    }
    if(lnum != readlen) {
        fprintf(stderr, "[E:%s] Unexpected index! Abort!\n", __func__);
        exit(EXIT_FAILURE);
    }
    free(buffer);
    return ret;
}

#endif // ARRAY_PARSER_H
