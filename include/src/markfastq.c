#include <zlib.h>
#include <stdio.h>
#include <regex.h>
#include <stdint.h>
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

#define KSEQ_TO_STRING(handle, read, index, pass) fprintf(handle, "%s ~#!#~|FP=%i|BS=%s\n%s\n+\n%s\n", read->name.s, pass, index->seq.s, read->seq.s, read->qual.s)


int main(int argc, char *argv[])
{
    FILE *outHandle;
    gzFile fp_read, fp_index;
    kseq_t *seq;
    kseq_t *seq_index;
    int l;
    int l_index;
    regex_t hp_regex;
    int regex_int_ret;
    int pass;
    char buffer [20];
    sprintf(buffer, "([ACTG])\1{%i}", 10);
    regex_int_ret = regcomp(&hp_regex, buffer, 0);
    if(regex_int_ret){
        fprintf(stderr, "Could not compile regex. Abort!\n");
        exit(1);
    }
    if (argc == 1 || argc < 3) {
        fprintf(stderr, "Usage: %s <Fq.seq> <Index.seq> <out.seq>"
                        "\nLeave out out.seq to emit to stdout", argv[0]);
        return 1;
    }
    if (argc == 3) {
        fprintf(stderr, "No path to output fastq provided. Defaulting to stdout.\n");
        outHandle = stdout;
    }
    else {
        outHandle = fopen(argv[3], "w");
    }
    fp_read = gzopen(argv[1], "r");
    fp_index = gzopen(argv[1], "r");
    seq = kseq_init(fp_read);
    seq_index = kseq_init(fp_index);
    while ((l = kseq_read(seq)) >= 0) {
        // Iterate through second fastq file.
        l_index = kseq_read(seq_index);
        if (l_index < 0) {
            fprintf(stderr, "Index return value for kseq_read less than 0. Are the fastqs different sizes? Abort!\n");
            return 1;
        }
        // Test regex.
        regex_int_ret = regexec(&hp_regex, seq_index->seq.s, 0, NULL, 0);
        // Set pass or fail for record.
        pass = regex_int_ret > 0 ? 0 : 1;
        KSEQ_TO_STRING(outHandle, seq, seq_index, pass);
    }
    fprintf("return value: %d\n", l);
    kseq_destroy(seq);
    gzclose(fp_read);
    fclose(outHandle);
    return 0;
}
