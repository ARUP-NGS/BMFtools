#include <zlib.h>
#include <stdio.h>
#include <stdint.h>
#include <regex.h>
#include <getopt.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include <stdlib.h>
#include <sys/resource.h>
#include "fqmarksplit.h"

// Inline function declarations
//int lh3_sort_call(char *fname, char *outfname);
void FREE_SPLITTER(mark_splitter_t var);
//void apply_lh3_sorts(sort_overlord_t *dispatcher, mss_settings_t *settings);
int ipow(int base, int exp);
mark_splitter_t init_splitter(mss_settings_t *settings_ptr);
int get_binner(char binner[], int length);
void splitmark_core(kseq_t *seq1, kseq_t *seq2, kseq_t *seq_index,
				    mss_settings_t settings, mark_splitter_t splitter);
sort_overlord_t build_mp_sorter(mark_splitter_t* splitter_ptr, mss_settings_t *settings_ptr);
void free_mp_sorter(sort_overlord_t var);
char test_hp(kseq_t *seq, int threshold);

// Macros



// Allocate file handle array memory, open file handles.

/*
#define FREE_SPLITTER(var) \
    for(int i_##var_tmp = 0; i_##var_tmp < var.n_handles; i_##var_tmp++) {\
        fclose(var.tmp_out_handles_r1[i_##var_tmp]);\
        fclose(var.tmp_out_handles_r2[i_##var_tmp]);\
        free(var.fnames_r1);\
        free(var.fnames_r2);\
    }\
    free(var.tmp_out_handles_r1);\
    free(var.tmp_out_handles_r2);\
    free(var.fnames_r1);\
    free(var.fnames_r2);
*/


void print_usage(char *argv[])
{
        fprintf(stderr, "Usage: %s <options> -i <Index.seq> <Fq.R1.seq> <Fq.R2.seq>"
                        "\nFlags:\n"
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

void print_opt_err(char *argv[], char *optarg)
{
    print_usage(argv);
    fprintf(stderr, "Unrecognized option %s. Abort!\n", optarg);
    exit(1);
}

int main(int argc, char *argv[])
{
    // Build settings struct
    int hp_threshold;
    int n_nucs;
    char *output_basename;
    int threads;
    const char *default_basename = "metasyntactic_var";
    mss_settings_t settings = {
        .hp_threshold = 10,
        .n_nucs = 4,
        .index_fq_path = NULL,
        .output_basename = NULL,
        .threads = 4,
        .input_r1_path = NULL,
        .input_r2_path = NULL,
        .n_handles = 0,
        .notification_interval = 100000
    };
    settings.n_handles = ipow(4, settings.n_nucs);
    int c;
    while ((c = getopt(argc, argv, "t:h:o:i:n:")) > -1) {
        switch(c) {
            case 'n': settings.n_nucs = atoi(optarg); break;
            case 't': settings.hp_threshold = atoi(optarg); break;
            case 'o': settings.output_basename = strdup(optarg);break;
            case 'i': settings.index_fq_path = strdup(optarg); break;
            case 'h': print_usage(argv); return 0;
            default: print_opt_err(argv, optarg);
        }
    }

    if(!settings.index_fq_path) {
        fprintf(stderr, "Index fastq required. See usage.\n");
        print_usage(argv);
        return 1;
    }

    if(!settings.output_basename) {
        fprintf(stderr, "Output basename required. See usage.\n");
        print_usage(argv);
        return 1;
    }
    fprintf(stderr, "About to get the read paths.\n");
    char r1_fq_buf[100];
    char r2_fq_buf[100];
    if(argc - 1 != optind + 1) {
        fprintf(stderr, "Both read 1 and read 2 fastqs are required. See usage.\n", argc, optind);
        print_usage(argv);
        return 1;
    }
    strcpy(r1_fq_buf, argv[optind]);
    strcpy(r2_fq_buf, argv[optind + 1]);
    mss_settings_t *settings_ptr = &settings;
    gzFile fp_read1, fp_read2, fp_index;
    fp_read1 = gzopen(r1_fq_buf, "r");
    fp_read2 = gzopen(r2_fq_buf, "r");
    fp_index = gzopen(settings.index_fq_path, "r");
    kseq_t *seq1;
    kseq_t *seq2;
    kseq_t *seq_index;
    int l1, l2, l_index;
    seq1 = kseq_init(fp_read1);
    seq2 = kseq_init(fp_read2);
    seq_index = kseq_init(fp_index);
    mark_splitter_t splitter = init_splitter(settings_ptr);
    mark_splitter_t *splitter_ptr = &splitter;
/*
    fprintf(stderr, "Hey, can I even read this fastq? %s, %s, %i", seq1->seq.s, seq1->qual.s, l);
    fprintf(stderr, "Hey, my basename is %s\n", settings.output_basename);
*/
    splitmark_core(seq1, seq2, seq_index,
    			   settings, splitter);
    //apply_lh3_sorts(&splitter, &settings);
    /*
    sort_overlord_t dispatcher = build_mp_sorter(splitter_ptr, settings_ptr);
    apply_lh3_sorts(&dispatcher, settings_ptr);
    free_mp_sorter(dispatcher);
    */
    FREE_SETTINGS(settings);
    FREE_SPLITTER(splitter);
    return 0;
}
