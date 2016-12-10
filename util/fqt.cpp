#include <zlib.h>
#include <inttypes.h>
#include <memory>

#include "htslib/kseq.h"
#include "htslib/kstring.h"
#include "dlib/logging_util.h"
KSEQ_INIT(gzFile, gzread)


int usage(char **argv, int rc) {
    fprintf(stderr, "Usage: %s -l<blen> infq1 infq2 outfq1 outfq2 outfq3\n For converting non-duplex inline barcodes for compatibility with collapse secondary", *argv);
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
    int blen(-1);
    int c;
    while((c = getopt(argc, argv, "l:h?")) > 0) {
        switch(c) {
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
    char obs[1024], obq[2014];
    obs[blen * 2] = obq[blen * 2] = 0;
    while(kseq_read(ks1) >= 0 && kseq_read(ks2) >= 0) {
        fqw(ks1, blen, out1);
        fqw(ks2, blen, out2);
        memcpy(obs, ks1->seq.s, blen);
        memcpy(obs + blen, ks2->seq.s, blen);
        memcpy(obq, ks1->qual.s, blen);
        memcpy(obq + blen, ks2->qual.s, blen);
        bcw(ks1, obs, obq, outindex);
    }
    kseq_destroy(ks1); kseq_destroy(ks2);
    gzclose(in1); gzclose(in2);
    fclose(out1); fclose(out2); fclose(outindex);
    return 0;
}
