#include "kseq.h"

KSEQ_INIT(gzFile, gzread)

// Typedefs
typedef struct markfastq_settings {
    int single_line;
    int hp_threshold;
    int write_to_stdout;
    int n_nucs;
    char * output_basename;
    char * index_fq_path;
} markfastq_settings_t;

typedef struct mark_splitter {
    FILE **out_handles;
    int n_nucs;
    int n_handles;
    char **fnames;
} mark_splitter_t;

int ipow(int base, int exp);
void write_kseq(FILE *handle, kseq_t *read, kseq_t *index, int pass, markfastq_settings_t *settings);
