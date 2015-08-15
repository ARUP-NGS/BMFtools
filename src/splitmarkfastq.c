#include <zlib.h>
#include <stdio.h>
#include <stdint.h>
#include <regex.h>
#include <getopt.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include "splitmarkfastq.h"
#include "include/sort/lh3sort.c"



// Pre-macro definitions


// Macros

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

#define SWITCH_FMODE_SPLITTER(var) for (int i = 0; i < var.n_handles; i++) \
    {\
        fclose(var.out_handles[i]);\
        var.out_handles[i] = fopen(var.fnames[i], "r");\
    }\


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
    size_t count = 0;\
    for(i_##binner = length_##binner;i_##binner > 0;i_##binner--){\
        char_to_num(binner[i_##binner - 1], inc_##binner);\
        bin += (ipow(4, count) * inc_##binner);\
        count++;\
    }

#define FREE_SETTINGS(settings) free(settings.output_basename);\
    free(settings.index_fq_path);

#define lh3_sort_call(fname, retvar)\
    int ret_lh3_sort_##fname;\
    char **lh3_argv_##fname = (char **)malloc(4 * sizeof(char *));\
    lh3_argv_##fname[1] = strdup("-t\'|\'");\
    lh3_argv_##fname[2] = strdup("-k2,2");\
    lh3_argv_##fname[3] = strdup(fname);\
    retvar = lh3_sort_main(4, lh3_argv_##fname);\
    free(lh3_argv_##fname[1]);\
    free(lh3_argv_##fname[2]);\
    free(lh3_argv_##fname[3]);\
    free(lh3_argv_##fname);

  
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
                        "\nFlags:\n"
                        "-f: Write each record as a single line. Default: True.\n"
                        "-t: Homopolymer failure threshold. A molecular barcode with"
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
    int l;
    int l_index;
    int pass;
    char buffer [20];
    char input_fastq_buffer [100];
    fprintf(stderr, "Hey, I started going in main.\n");

    // Build settings struct
    markfastq_settings_t settings = {
        .single_line = 1,
        .hp_threshold = 10,
        .n_nucs = 1,
        .index_fq_path = NULL,
        .output_basename = NULL,
        .threads = 1
    };
    // Parse in command-line options.
    int c;
    fprintf(stderr, "Hey, I started parsing my options.\n");
    while ((c = getopt(argc, argv, "fh:o:i:n:p:")) > -1) {
        switch(c) {
            case 'f': settings.single_line = 0; break;
            case 'n': settings.n_nucs = atoi(optarg); break;
            case 't': settings.hp_threshold = atoi(optarg); break;
            case 'o': settings.output_basename = strdup(optarg);break;
            case 'i': settings.index_fq_path = strdup(optarg); break;
            case 'p': settings.threads = atoi(optarg); break;
            case 'h': print_usage(argv); return 0;
            default: print_opt_err(argv);
        }
    }
    fprintf(stderr, "Hey, my basename is %s\n", settings.output_basename);

    if(!settings.index_fq_path) {
        fprintf(stderr, "Index fastq required. See usage.\n");
        print_usage(argv);
        return 1;
    }

    if(!settings.index_fq_path) {
        fprintf(stderr, "Index fastq required. See usage.\n");
        print_usage(argv);
        return 1;
    }

    //fprintf(stderr, "Hey, I parsed my options.\n");

    int numHandles = ipow(4, settings.n_nucs);
    //fprintf(stderr, "Hey, i have %i nucs set.", settings.n_nucs);

    markfastq_settings_t *settings_p = &settings;
    // Build file handle struct
    INIT_SPLITTER(splitter, settings_p)
    fprintf(stderr, "Hey, I initialized my file handles.\n");
    // Build regex
    sprintf(buffer, "([ACTG])\1{%i,}", settings.hp_threshold);
    regex_t regex;
    int reti;

    reti = regcomp(&regex, buffer, REG_EXTENDED);
    if(reti) {
        fprintf(stderr, "Could not compile regular expression '%s'.\n", buffer);
        exit(1);
    }
    strcpy(input_fastq_buffer, argv[optind]);
    fp_read = gzopen(input_fastq_buffer, "r");
    fp_index = gzopen(settings.index_fq_path, "r");
    seq = kseq_init(fp_read);
    seq_index = kseq_init(fp_index);
    char binner [10];
    int bin_num;
    while ((l = kseq_read(seq)) >= 0) {
        // Iterate through second fastq file.
        l_index = kseq_read(seq_index);
        if (l_index < 0) {
            fprintf(stderr, "Index return value for kseq_read less than "
                            "0. Are the fastqs different sizes? Abort!\n");
            return 1;
        }
        // Set pass or fail for record.
        reti = regexec(&regex, seq_index->seq.s, 0, NULL, 0);
        pass = (reti) ? 1 : 0;
        if(!pass){
            printf("Is this right? Should string %s have been failed?\n", seq_index->qual);
            return 1;
        }
        strncpy(binner, seq_index->seq.s, settings.n_nucs);
        binner[settings.n_nucs] = '\0';
        FIND_BIN(binner, bin_num)
        write_kseq(splitter.out_handles[bin_num],
                   seq, seq_index, pass, &settings);
    }
    kseq_destroy(seq);
    gzclose(fp_read);
    FREE_SPLITTER(splitter)
    FREE_SETTINGS(settings)
    return 0;
}
