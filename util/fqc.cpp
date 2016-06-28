#include <zlib.h>
#include <inttypes.h>

#include "htslib/kseq.h"
#include "htslib/kstring.h"
#include "dlib/logging_util.h"
KSEQ_INIT(gzFile, gzread)

int main(int argc, char **argv)
{
    if(argc != 2) LOG_EXIT("Usage: %s <inpath>\n", argv[0]);
    gzFile in(gzopen(argv[1], "r"));
    if(in == nullptr) LOG_EXIT("Could not open file %s.\n", argv[1]);
    kseq_t *seq(kseq_init(in));
    if(seq == nullptr) LOG_EXIT("Could not open file %s as fastq.\n", argv[2]);
    uint64_t c(0);
    for(;kseq_read(seq) >= 0; ++c);
    fprintf(stderr, "#Filename\tCount\n%s\t%lu\n", argv[1], c);
    kseq_destroy(seq);
    gzclose(in);
}
