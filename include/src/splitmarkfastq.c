#include <zlib.h>
#include <stdio.h>
#include <stdint.h>
#include <regex.h>
#include <getopt.h>
#include <math.h>
#include <string.h>
#include "kseq.h"

// Macros

// Print fastq record in standard format. (4 lines per record)
#define KSEQ_TO_STRING(handle, read, index, pass) fprintf(handle, "%s ~#!#~|FP=%i|BS=%s\n%s\n+\n%s\n", read->name.s, pass, index->seq.s, read->seq.s, read->qual.s)
// Print fastq record in single line format. (1 line per record, fields separated by tabs. Used for use cases involving GNU sort.)
#define KSEQ_TO_SINGLE_LINE(handle, read, index, pass) fprintf(handle, "%s ~#!#~|FP=%i|BS=%s\t%s\t+\t%s\n", read->name.s, pass, index->seq.s, read->seq.s, read->qual.s)

// Allocate file handle array memory, open file handles.
#define INIT_SPLITTER(var, n_nucs, basename) mark_splitter_t var = {\
    .n_nucs = n_nucs;\
    .out_handles = (FILE **)malloc(ipow(4, n_nucs) * sizeof(FILE *));\
    .n_handles = ipow(4, n_nucs);\
    .fnames = (char **)malloc(ipow(4, n_nucs) * sizeof(char *));\
}\
char tmp_##var##_buffer [100];\
for (int i = 0; i < var.n_handles; i++) {\
    sprintf(tmp_##var##_buffer, "%s.%i.fastq", basename, i);\
    size_t length = strlen(tmp_##var##_buffer);\
    var.fnames[i] = (char *)malloc(i(strlen(tmp_##var##buffer) + 1) * sizeof(char));\
    strcpy(var.fnames[i], tmp_##var##_buffer);\
    var.out_handles[i] = fopen(var.fnames[i], "w");\
}


// Close file handles, then free the malloc'd array.
#define FREE_SPLITTER(var) for(int i = 0; i < var.n_handles; i++) \
    {\
        fclose(var.out_handles[i]);\
        free(var.fnames[i]);\
    }\
    free(var.out_handles);

// Typedefs
typedef struct markfastq_settings {
    int single_line;
    int hp_threshold;
    int write_to_stdout;
    int number_nucs;
    char * output_basename;
    char * index_fq_path;
} markfastq_settings_t;

typedef struct mark_splitter {
    FILE **out_handles;
    int n_nucs;
    int n_handles;
    char **fnames;
} mark_splitter_t;
  
// Functions

inline void write_kseq(FILE *handle, kseq_t *read, kseq_t *index, int pass, markfastq_settings_t *settings)
{
    if(settings->single_line) {
        KSEQ_TO_SINGLE_LINE(handle, read, index, pass);
    }
    else {
        KSEQ_TO_STRING(handle, read, index, pass);
    }
}

inline int ipow(int base, int exp)
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

void print_usage()
{
        fprintf(stderr, "Usage: %s <options> -i <Index.seq> <Fq.seq>"
                        "\nFlags:\n-s: Write to stdout. Conflicts with n != 0.\n"
                        "-f: Write each record as a single line. Default: True.\n"
                        "-h: Homopolymer failure threshold. A molecular barcode with"
                        " a homopolymer of length >= this limit is flagged as QC fail.\n"
                        "-o: Output basename. Currently required, as string "
                        "manipulation in C is a bit of work and I'd rather spend my "
                        "time building code than messing around with string "
                        "manipulation that doesn't add to the code base."
                        "\n-i: Index fastq path. Required.\n"
                        "-n: Number of nucleotides at the beginning of the barcode to use to split the output.\n", argv[0]);
}

void print_opt_err()
{
    print_usage();
    fprintf(stderr, "Unrecognized option. Abort!\n");
    return 1;
}

KSEQ_INIT(gzFile, gzread)

int main(int argc, char *argv[])
{
    gzFile fp_read, fp_index;
    kseq_t *seq;
    kseq_t *seq_index;
    char * pch;
    int l;
    int l_index;
    int pass;
    char buffer [20];
    char basename_buffer [100];

    // Build settings struct
    markfastq_settings_t settings = {
        .single_line = 1;
        .hp_threshold = 10;
        .write_to_stdout = false;
        .number_nucs = 0;
        .index_fq_path = NULL;
        .output_basename = NULL;
    }
    // Parse in command-line options.
    while ((c = getopt(argc, argv, "sfnph:")) > -1) {
        switch(c) {
            case 's': settings.write_to_stdout = 1; break;
            case 'f': settings.single_line = 0; break;
            case 'n': settings.number_nucs = atoi(optarg); break;
            case 'h': settings.hp_threshold = atoi(optarg); break;
            case 'o': settings.output_basename = strcpy(optarg); break;
            case 'i': settings.index_fq_path = strcpy(optarg); break;
            default: print_opt_err();
        }
    }
    if(!settings.index_fq_path) {
        fprintf(stderr, "Index fastq required. See usage.\n");
        print_usage();
    }

    int numHandles = ipow(4, settings.number_nucs)

    // Build file handle struct
    INIT_SPLITTER(splitter, settings.n_nucs, settings.output_basename)

    // Build regex
    sprintf(buffer, "([ACTG])\1{%i,}", settings.hp_threshold);
    regex_t regex;
    int reti;

    reti = regcomp(&regex, buffer, REG_EXTENDED);
    if(reti) {
        fprintf(stderr, "Could not compile regular expression '%s'.\n", buffer);
        exit(1);
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
    char * binner = malloc((settings.n_nucs + 1) * sizeof(char));
    while ((l = kseq_read(seq)) >= 0) {
        // Iterate through second fastq file.
        l_index = kseq_read(seq_index);
        if (l_index < 0) {
            fprintf(stderr, "Index return value for kseq_read less than 0. Are the fastqs different sizes? Abort!\n");
            return 1;
        }
        // Set pass or fail for record.
        reti = regexec(&regex, seq_index->seq.s, 0, NULL, 0);
        pass = (reti) ? 1 : 0;
        if(!pass){
            printf("Is this right? Should string %s have been failed?\n", seq_index->qual);
            return 1;
        }
        memcpy(binner, &seq_index->seq.s, 2);
        write_kseq(splitter.out_handles[FIND_BIN(binner)],
                   seq, seq_index, pass, &settings);
    }
    fprintf("return value: %d\n", l);
    kseq_destroy(seq);
    gzclose(fp_read);
    fclose(outHandle);
    FREE_SPLITTER(splitter);
    return 0;
}
