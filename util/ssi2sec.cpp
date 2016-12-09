#include <zlib.h>
#include <inttypes.h>
#include <memory>

#include "htslib/kseq.h"
#include "htslib/kstring.h"
#include "dlib/logging_util.h"
KSEQ_INIT(gzFile, gzread)


int usage(char **argv, int rc) {
    fprintf(stderr, "Usage: %s -l<blen> infq1 infq2 outfq1 outfq2 outfq3\n For converting non-duplex inline barcodes for compatibility with collapse secondary\n"
                    "-2: Get single-sided barcode from read2.\n", *argv);
    return rc;
}


void fqw(kseq_t *ks, int offset, FILE *fp) {
    fprintf(fp, "@%s %s\n%s\n+\n%s\n", ks->name.s,
            ks->comment.s, ks->seq.s + offset, ks->qual.s + offset);
}

void bcw(kseq_t *ks, char *s, char *q, FILE *fp) {
    fprintf(fp, "@%s %s\n%s\n+\n%s\n", ks->name.s,
            ks->comment.s, s, q);
}


int main(int argc, char **argv)
{
    if(argc < 6) return usage(argv, EXIT_FAILURE);
    int blen(-1), use_read1(1);
    int c;
    while((c = getopt(argc, argv, "l:2h?")) > 0) {
        switch(c) {
        case '2': use_read1 = 0;
        case 'l': blen  = atoi(optarg); break;
        case 'h': case '?': return usage(argv, EXIT_SUCCESS);
        }
    }
    if(blen < 0) {
        fprintf(stderr, "blen required. Abort!\n");
        return usage(argv, EXIT_FAILURE);
    }
    gzFile in1(gzopen(argv[optind], "rb")), in2(gzopen(argv[optind + 1], "rb"));
    kseq_t *ks1(kseq_init(in1)), *ks2(kseq_init(in2));
    FILE *out1(fopen(argv[optind + 2], "w")), *out2(fopen(argv[optind + 3], "w")), *outindex(fopen(argv[optind + 4], "w"));
    char bpos[128];
    bpos[blen] = '\0';
    kseq_t *bks(use_read1 ? ks1: ks2);
    while(kseq_read(ks1) >= 0 && kseq_read(ks2) >= 0) {
        memcpy(bpos, bks->seq.s, blen);
        fprintf(outindex, "%s barcode\n%s\n+\n", ks1->name.s, bpos);
        fwrite(bks->qual.s, 20, 1, outindex);
        fputc('\n', outindex);
        memset(bks->seq.s, 'N', blen);
        memset(bks->qual.s, '#', blen);
        fqw(ks1, 0, out1);
        fqw(ks1, 0, out2);
    }
    kseq_destroy(ks1); kseq_destroy(ks2);
    gzclose(in1); gzclose(in2);
    fclose(out1); fclose(out2); fclose(outindex);
    return 0;
}
