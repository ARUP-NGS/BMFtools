#include "include/kseq.h"

KSEQ_INIT(gzFile, gzread)

// Typedefs
typedef struct markfastq_settings {
    int single_line;
    int hp_threshold;
    int n_nucs;
    char * output_basename;
    int threads;
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

// Print fastq record in standard format. (4 lines per record)                                                   
#define KSEQ_TO_STRING(handle, read, index, pass) fprintf(handle,\
        "%s FP:i:%i|BS:Z:%s\n%s\n+\n%s\n",\
        read->name.s, pass, index->seq.s, read->seq.s, read->qual.s)
// Print fastq record in single line format. (1 line per record, fields separated by tabs. Used for use cases involving GNU sort.)
#define KSEQ_TO_SINGLE_LINE(handle, read, index, pass) fprintf(handle,\
        "%s FP:i:%i|BS:Z:%s\t%s\t+\t%s\n",\
        read->name.s, pass, index->seq.s, read->seq.s, read->qual.s)
