#include <zlib.h>
#include <stdio.h>
#include <stdint.h>
#include <regex.h>
#include <getopt.h>
#include <math.h>
#include <string.h>
#include "kseq.h"
#include "splitmarkfastq.h"



// Pre-macro definitions


// Macros
// Print fastq record in standard format. (4 lines per record)
#define KSEQ_TO_STRING(handle, read, index, pass) fprintf(handle, "%s ~#!#~|FP=%i|BS=%s\n%s\n+\n%s\n", read->name.s, pass, index->seq.s, read->seq.s, read->qual.s)
// Print fastq record in single line format. (1 line per record, fields separated by tabs. Used for use cases involving GNU sort.)
#define KSEQ_TO_SINGLE_LINE(handle, read, index, pass) fprintf(handle, "%s ~#!#~|FP=%i|BS=%s\t%s\t+\t%s\n", read->name.s, pass, index->seq.s, read->seq.s, read->qual.s)

// Allocate file handle array memory, open file handles.
#define INIT_SPLITTER(var, settings_ptr) mark_splitter_t var = {\
    .out_handles = (FILE **)malloc(ipow(4, settings_ptr->n_nucs) * sizeof(FILE *)),\
    .n_handles = ipow(4, settings_ptr->n_nucs),\
    .fnames = (char **)malloc(ipow(4, settings_ptr->n_nucs) * sizeof(char *)),\
    .n_nucs = settings_ptr->n_nucs,\
};\
char tmp_##var##_buffer [100];\
for (int i = 0; i < var.n_handles; i++) {\
    sprintf(tmp_##var##_buffer, "%s.%i.fastq", settings_ptr->output_basename, i);\
    size_t length = strlen(tmp_##var##_buffer);\
    var.fnames[i] = (char *)malloc((strlen(tmp_##var##_buffer) + 1) * sizeof(char));\
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
// Select which FILE ** index should be used to write this read from a char *.

#define char_to_num(character, increment) switch(character) {\
    case 'C' : increment = 1; break;\
    case 'G' : increment = 2; break;\
    case 'T' : increment = 3; break;\
    default: increment = 0; break;\
    }

#define FIND_BIN(binner, bin) bin = 0;\
    size_t length_##binner = strlen(binner);\
    int inc_##binner;\
    int i_##binner = length_##binner;\
    char tmpchar;\
    while (i_##binner){\
        char_to_num(binner[--i_##binner], inc_##binner);\
        bin += inc_##binner;\
    }

  
// Functions

void write_kseq(FILE *handle, kseq_t *read, kseq_t *index, int pass, markfastq_settings_t *settings)
{
    if(settings->single_line) {
        KSEQ_TO_SINGLE_LINE(handle, read, index, pass);
    }
    else {
        KSEQ_TO_STRING(handle, read, index, pass);
    }
}

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

void print_usage(char *argv[])
{
        fprintf(stderr, "Usage: %s <options> -i <Index.seq> <Fq.seq>"
                        "\nFlags:\n-s: Write to stdout. Conflicts with n != 0.\n"
                        "-f: Write each record as a single line. Default: True.\n"
                        "-h: Homopolymer failure threshold. A molecular barcode with"
                        " a homopolymer of length >= this limit is flagged as QC fail."
                        "Default: 10.\n"
                        "-o: Output basename. Currently required, as string "
                        "manipulation in C is a bit of work and I'd rather spend my "
                        "time building code than messing around with string "
                        "manipulation that doesn't add to the code base.\n"
                        "-i: Index fastq path. Required.\n"
                        "-n: Number of nucleotides at the beginning of the barcode to use to split the output.\n", argv[0]);
}

void print_opt_err(char *argv[])
{
    print_usage(argv);
    fprintf(stderr, "Unrecognized option. Abort!\n");
    exit(1);
}


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
        .single_line = 1,
        .hp_threshold = 10,
        .write_to_stdout = 0,
        .n_nucs = 0,
        .index_fq_path = NULL,
        .output_basename = NULL,
    };
    // Parse in command-line options.
    int c;
    while ((c = getopt(argc, argv, "sfnph:")) > -1) {
        switch(c) {
            case 's': settings.write_to_stdout = 1; break;
            case 'f': settings.single_line = 0; break;
            case 'n': settings.n_nucs = atoi(optarg); break;
            case 'h': settings.hp_threshold = atoi(optarg); break;
            case 'o': strcpy(settings.output_basename, optarg); break;
            case 'i': strcpy(settings.index_fq_path, optarg); break;
            default: print_opt_err(argv);
        }
    }
    if(!settings.index_fq_path) {
        fprintf(stderr, "Index fastq required. See usage.\n");
        print_usage(argv);
    }

    int numHandles = ipow(4, settings.n_nucs);

    markfastq_settings_t *settings_p = &settings;
    // Build file handle struct
    INIT_SPLITTER(splitter, settings_p)

    // Build regex
    sprintf(buffer, "([ACTG])\1{%i,}", settings.hp_threshold);
    regex_t regex;
    int reti;

    reti = regcomp(&regex, buffer, REG_EXTENDED);
    if(reti) {
        fprintf(stderr, "Could not compile regular expression '%s'.\n", buffer);
        exit(1);
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
        int bin_num;
        FIND_BIN(binner, bin_num)
        write_kseq(splitter.out_handles[bin_num],
                   seq, seq_index, pass, &settings);
    }
    kseq_destroy(seq);
    gzclose(fp_read);
    FREE_SPLITTER(splitter);
    return 0;
}
