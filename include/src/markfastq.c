#include <zlib.h>
#include <stdio.h>
#include <stdint.h>
#include <regex.h>
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

#define KSEQ_TO_STRING(handle, read, index, pass) fprintf(handle, "%s ~#!#~|FP=%i|BS=%s\n%s\n+\n%s\n", read->name.s, pass, index->seq.s, read->seq.s, read->qual.s)

#define test_hp(read, regex)

int main(int argc, char *argv[])
{
    FILE *outHandle;
    gzFile fp_read, fp_index;
    kseq_t *seq;
    kseq_t *seq_index;
    char * pch;
    int l;
    int l_index;
    int pass;
    char buffer [20];
    const int hp_threshold = 10;
    sprintf(buffer, "([ACTG])\1{%i,}", hp_threshold);
    printf("Please run.\n");

    // Build regex
    regex_t regex;
    int reti;

    reti = regcomp(&regex, buffer, 0);
    if(reti) {
        fprintf(stderr, "Could not compile regular expression '%s'.\n", buffer);
        exit(1);
    }

    if (argc == 1 || argc < 3) {
        fprintf(stderr, "Usage: %s <Fq.seq> <Index.seq> <out.seq>"
                        "\nLeave out out.seq to emit to stdout.\n", argv[0]);
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
    fp_index = gzopen(argv[2], "r");
    seq = kseq_init(fp_read);
    seq_index = kseq_init(fp_index);
    while ((l = kseq_read(seq)) >= 0) {
        // Iterate through second fastq file.
        l_index = kseq_read(seq_index);
        if (l_index < 0) {
            fprintf(stderr, "Index return value for kseq_read less than 0. Are the fastqs different sizes? Abort!\n");
            return 1;
        }
        // Set pass or fail for record.
        reti = regexec(&regex, seq_index->qual.s, 0, NULL, 0);
        pass = (reti) ? 1 : 0;
        if(!pass){
            printf("Is this right? Should string %s have been failed?\n", seq_index->qual);
            exit(1);
        }
        KSEQ_TO_STRING(outHandle, seq, seq_index, pass);
    }
    fprintf("return value: %d\n", l);
    kseq_destroy(seq);
    gzclose(fp_read);
    fclose(outHandle);
    return 0;
}
