#include "include/kseq.h"

KSEQ_INIT(gzFile, gzread)

// Typedefs
typedef struct mss_settings {
    int hp_threshold;
    int n_nucs;
    char *output_basename;
    int threads;
    char *index_fq_path;
} mss_settings_t;

typedef struct mark_splitter {
    FILE **tmp_out_handles;
    int n_nucs;
    int n_handles;
    char **fnames;
} mark_splitter_t;

typedef struct sort_overlord {
    mark_splitter_t splitter;
    FILE **sort_out_handles;
    char **out_fnames;
} sort_overlord_t;

int ipow(int base, int exp)
{
    int result = 1;
    while (exp)
    {
        if (exp & 1)
            result *= base;
        exp >>= 1;
        base *= base;
    }

    return result;
}

int lh3_sort_call(char *fname, char *outfname);
void FREE_SPLITTER(mark_splitter_t var);
void apply_lh3_sorts(sort_overlord_t *dispatcher, mss_settings_t *settings);

// Print fastq record in single line format. (1 line per record, fields separated by tabs. Used for use cases involving GNU sort.)
#define KSEQ_TO_SINGLE_LINE(handle, read, index, pass) fprintf(handle,\
        "%s FP:i:%i|BS:Z:%s\t%s\t+\t%s\n",\
        read->name.s, pass, index->seq.s, read->seq.s, read->qual.s)
