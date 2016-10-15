#ifndef ARRAY_PARSER_H
#define ARRAY_PARSER_H
#include <assert.h>
#include "dlib/io_util.h"

#define NQSCORES 45uL // Number of q scores in sequencing.

#ifdef __GNUC__
#define INLINE __attribute__((always_inline)) inline
#else
#define INLINE inline
#endif

namespace bmf {

#if !NDEBUG
// The rescale_qscore is a little obfuscated. A human-readable equivalent is included before it.
// In debug mode, it asserts that both versions return the same value.
CONST static INLINE char rescale_qscore_check(int readnum, char qscore, int cycle, char base, int readlen, char *rescaler)
{
    if(base == 'N') return '#';
    int index(readnum);
    int mult(2);
    index += cycle * mult;
    mult *= readlen;
    index += (qscore - 35) * mult; // Subtract 35 - 33 to get to phred space, 2 to offset by 2.
    mult *= NQSCORES;
    index += mult * nuc2num_acgt(base);
    return rescaler[index] + 33;
}
#endif // !NDEBUG for rescale_qscore_check.

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
CONST static INLINE char rescale_qscore(int readnum, char qscore, int cycle, char base, int readlen, char *rescaler) {
    return (base == 'N') ? '#'
                         : rescaler[readnum + (cycle << 1) + (readlen << 1) * (qscore - 35) + (NQSCORES * (readlen << 1)) * nuc2num_acgt(base)] + 33;
}

static void period_to_null(char *c)
{
    for(;*c;++c) if(*c == '.') *c = 0;
}


static inline char *parse_1d_rescaler(char *qual_rescale_fname)
{
    int length, lnum;
    size_t arr_len, index;
    FILE *fp;
    char *buffer, *ret, *tok;
    const int readlen(dlib::count_lines(qual_rescale_fname));
    LOG_DEBUG("Number of lines: %i.\n", readlen);
    if((fp = fopen(qual_rescale_fname, "rb")) == nullptr) {
        LOG_EXIT("Could not open file %s. Abort mission!\n", qual_rescale_fname);
    }
    // Get the length of the file contents.
    fseek(fp, 0, SEEK_END);
    length = ftell(fp);
    fseek(fp, 0, SEEK_SET);
    buffer = (char *)malloc((size_t)length);
    if(!buffer) LOG_EXIT("Could not allocate memory.\n");
    fread(buffer, 1, length, fp);
    fclose(fp), fp = nullptr;
    arr_len = 2 * readlen * NQSCORES * 4;
    ret = (char *)malloc(arr_len * sizeof(char));
    memset(ret, -127, arr_len); // Set all of these char values to -127, which is definitely unprintable
    tok = nullptr;
    index = 0;
    LOG_DEBUG("Parsing in array with read len %lu from %s...\n",
              arr_len, qual_rescale_fname);
    for(lnum = 0; lnum < readlen; ++lnum) {
        for(int readnum = 0; readnum < 2; ++readnum) {
            for(unsigned qnum = 0; qnum < NQSCORES; ++qnum) {
                for(int bnum = 0; bnum < 4; ++bnum) {
                    if((tok = strtok(tok ? nullptr : buffer, "|:,\n")) == nullptr)
                        break;
                    period_to_null(tok);
                    ret[index++] = static_cast<char>(atof(tok) + 0.5); // Round up.
                    if(ret[index - 1] > 93) {
                        LOG_WARNING("Rescaled quality score above the max"
                                    " that can be held in the fastq format. Capping at 93. Value: %i.\n",
                                    ret[index - 1]);
                        ret[index - 1] = 93;
                    } else if(ret[index - 1] < 2) ret[index - 1] = 2;
                }
            }
        }
    }
    assert(index == arr_len);
    assert(lnum == readlen);
    free(buffer);
    return ret;
}

} /* namespace bmf */

#endif // ARRAY_PARSER_H
