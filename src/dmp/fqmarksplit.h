#include <zlib.h>
#include <stdio.h>
#include <stdint.h>
#include <regex.h>
#include <getopt.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include "include/kseq.h"

#ifndef METASYNTACTIC_FNAME_BUFLEN
#define METASYNTACTIC_FNAME_BUFLEN 100
#endif

KSEQ_INIT(gzFile, gzread)


typedef struct mss_settings {
    int hp_threshold;
    int n_nucs;
    int n_handles;
    char *output_basename;
    int threads;
    int notification_interval;
    char *index_fq_path;
    char *input_r1_path;
    char *input_r2_path;
} mss_settings_t;

typedef struct mark_splitter {
    FILE **tmp_out_handles_r1;
    FILE **tmp_out_handles_r2;
    int n_nucs;
    int n_handles;
    char **fnames_r1;
    char **fnames_r2;
} mark_splitter_t;


typedef struct sort_overlord {
    mark_splitter_t *splitter;
    //FILE **sort_out_handles_r1;
    //FILE **sort_out_handles_r2;
    char **out_fnames_r1;
    char **out_fnames_r2;
} sort_overlord_t;

int ipow(int base, int exp);
void FREE_SPLITTER(mark_splitter_t var);
//void apply_lh3_sorts(sort_overlord_t *dispatcher, mss_settings_t *settings);

#ifndef KSEQ_2_FQ
#define KSEQ_2_FQ(handle, read, index, pass) fprintf(handle,\
        "@%s ~#!#~|FP=%i|BS=%s\n%s\n+\n%s\n",\
    read->name.s, pass, index->seq.s, read->seq.s, read->qual.s)
#endif


#ifndef FREE_SETTINGS
#define FREE_SETTINGS(settings) free(settings.output_basename);\
    free(settings.index_fq_path);\
    free(settings.input_r1_path);\
    free(settings.input_r2_path);
#endif

inline int test_hp(kseq_t *seq, int threshold)
{
	int run = 0;
	char last = '\0';
	for(int i = 0; i < seq->seq.l; i++){
		if(seq->seq.s[i] == 'N') {
			return 0;
		}
		if(seq->seq.s[i] == last) {
			run += 1;
		}
		else {
			run = 0;
			last = seq->seq.s[i];
		}
	}
	return (run < threshold);
}
/*
inline sort_overlord_t build_mp_sorter(mark_splitter_t* splitter_ptr, mss_settings_t *settings_ptr)
{
    sort_overlord_t ret = {
            .splitter = splitter_ptr,
            //.sort_out_handles_r1 = (FILE **)malloc(splitter_ptr->n_handles * sizeof(FILE *)),
            //.sort_out_handles_r2 = (FILE **)malloc(splitter_ptr->n_handles * sizeof(FILE *)),
            .out_fnames_r1 = (char **)malloc(ipow(4, settings_ptr->n_nucs) * sizeof(char *)),
            .out_fnames_r2 = (char **)malloc(ipow(4, settings_ptr->n_nucs) * sizeof(char *))
        };
    char tmp_buffer [METASYNTACTIC_FNAME_BUFLEN];
    size_t length;
    for (int i = 0; i < splitter_ptr->n_handles; i++) {
        sprintf(tmp_buffer, "%s.R1.%i.sort.fastq", settings_ptr->output_basename, i);
        length = strlen(tmp_buffer);
        ret.out_fnames_r1[i] = strdup(tmp_buffer);
        //ret.sort_out_handles_r1[i] = fopen(ret.out_fnames_r1[i], "w");
        sprintf(tmp_buffer, "%s.R2.%i.sort.fastq", settings_ptr->output_basename, i);
        ret.out_fnames_r2[i] = strdup(tmp_buffer);
        //ret.sort_out_handles_r2[i] = fopen(ret.out_fnames_r2[i], "w");
    }
    return ret;
}

inline void free_mp_sorter(sort_overlord_t var){
    for (int i = 0; i < var.splitter->n_handles; i++) {
        //fclose(var.sort_out_handles_r1[i]);
        //fclose(var.sort_out_handles_r2[i]);
        free(var.out_fnames_r1[i]);
        free(var.out_fnames_r2[i]);
    }
    //free(var.sort_out_handles_r1);
    //free(var.sort_out_handles_r2);
    free(var.out_fnames_r1);
    free(var.out_fnames_r2);
    //FREE_SPLITTER(*var.splitter);
    return;
}
*/


// Functions
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
/*
inline int lh3_sort_call(char *fname, char *outfname)
{
	char buffer[200];
    int retvar;
    char **lh3_argv = (char **)malloc(6 * sizeof(char *));
    lh3_argv[0] = strdup("lh3sort");
    lh3_argv[1] = strdup("-t'|'");
    lh3_argv[2] = strdup("-k2,2");
    lh3_argv[3] = strdup("-o");
    lh3_argv[4] = strdup(outfname);
    lh3_argv[5] = strdup(fname);
    sprintf(buffer, "%s %s %s %s %s %s", lh3_argv[0], lh3_argv[1], lh3_argv[2], lh3_argv[3], lh3_argv[4], lh3_argv[5]);
    fprintf(stderr, buffer + '\n');
    system(buffer);
    //retvar = lh3_sort_main(6, lh3_argv);
    for(int i = 1; i < 6; i++) {
        free(lh3_argv[i]);
    }
    free(lh3_argv);
    return retvar;
}
*/


inline void FREE_SPLITTER(mark_splitter_t var){
    for(int i = 0; i < var.n_handles; i++)
    {
        fclose(var.tmp_out_handles_r1[i]);
        fclose(var.tmp_out_handles_r2[i]);
        free(var.fnames_r1[i]);
        free(var.fnames_r2[i]);
    }
    free(var.tmp_out_handles_r1);
    free(var.tmp_out_handles_r2);
    return;
}
/*
void apply_lh3_sorts(sort_overlord_t *dispatcher, mss_settings_t *settings)
{
    int abort = 0;
    int index = -1;
    omp_set_num_threads(settings->threads);
    fprintf(stderr, "Number of threads: %i.\n", settings->threads);
    fprintf(stderr, "Number of : handles %i.\n", dispatcher->splitter->n_handles * 4);
    #pragma omp parallel for
    for(int i = 0; i < dispatcher->splitter->n_handles; i++) {
    	fprintf(stderr, "Now about to call an lh3 sort # %i. Input: %s. Output: %s.\n", i, dispatcher->splitter->fnames_r1[i], dispatcher->out_fnames_r1[i]);
        #pragma omp flush(abort)
		fprintf(stderr, "About to try opening file %s.\n", dispatcher->out_fnames_r1[i]);
        int ret = lh3_sort_call(dispatcher->splitter->fnames_r1[i], dispatcher->out_fnames_r1[i]);
        fprintf(stderr, "lh3_sort_call return value: %i.\n", ret);
        if(ret) {
            abort = 1;
            index = i;
            #pragma omp flush (abort)
            #pragma omp flush (index)
        }
    }
    if(abort) {
        fprintf(stderr,
                "lh3 sort call failed for file handle %s. (Non-zero exit status). Abort!",
                dispatcher->splitter->fnames_r1[index]);
                //free_mp_sorter(*dispatcher); // Delete allocated memory.
                // Will need to rewrite this for paired-end.
        exit(EXIT_FAILURE);
    }
    abort = 0;
    index = -1;
    #pragma omp parallel for
    for(int i = 0; i < dispatcher->splitter->n_handles; i++) {
        #pragma omp flush(abort)
    	fprintf(stderr, "Now about to call an lh3 sort # %i. Input: %s. Output: %s.\n", i, dispatcher->splitter->fnames_r2[i], dispatcher->out_fnames_r2[i]);
        int ret = lh3_sort_call(dispatcher->splitter->fnames_r2[i], dispatcher->out_fnames_r2[i]);
        if(ret) {
            abort = 1;
            index = i;
            #pragma omp flush (abort)
            #pragma omp flush (index)
        }
    }
    if(abort) {
        fprintf(stderr,
                "lh3 sort call failed for file handle %s. (Non-zero exit status). Abort!",
                dispatcher->splitter->fnames_r2[index]);
                //free_mp_sorter(*dispatcher); // Delete allocated memory.
                // Will need to rewrite this for paired-end.
        exit(EXIT_FAILURE);
    }
    return;
}
*/

#define char_to_num(character, increment) switch(character) {\
    case 'C' : increment = 1; break;\
    case 'G' : increment = 2; break;\
    case 'T' : increment = 3; break;\
    default: increment = 0; break;\
    }


inline int get_binner(char *barcode, int length) {
    int bin = 0;
    int inc_binner;
    size_t count = 0;
    for(int i = length;i > 0;i--){
        char_to_num(barcode[i - 1], inc_binner);
        bin += (ipow(4, count) * inc_binner);
        count++;
    }
    return bin;
}

inline mark_splitter_t init_splitter(mss_settings_t* settings_ptr)
{
    mark_splitter_t ret = {
        .n_handles = ipow(4, settings_ptr->n_nucs),
        .n_nucs = settings_ptr->n_nucs,
        .fnames_r1 = NULL,
        .fnames_r2 = NULL,
        .tmp_out_handles_r1 = NULL,
        .tmp_out_handles_r2 = NULL
    };
    // Avoid passing more variables.
    ret.tmp_out_handles_r1 = (FILE **)malloc(settings_ptr->n_handles * sizeof(FILE *));
    ret.tmp_out_handles_r2 = (FILE **)malloc(settings_ptr->n_handles * sizeof(FILE *));
    ret.fnames_r1 = (char **)malloc(ret.n_handles * sizeof(char *));
    ret.fnames_r2 = (char **)malloc(ret.n_handles * sizeof(char *));
    char tmp_buffer [METASYNTACTIC_FNAME_BUFLEN];
    for (int i = 0; i < ret.n_handles; i++) {
        sprintf(tmp_buffer, "%s.tmp.%i.R1.fastq", settings_ptr->output_basename, i);
        ret.fnames_r1[i] = strdup(tmp_buffer);
        sprintf(tmp_buffer, "%s.tmp.%i.R2.fastq", settings_ptr->output_basename, i);
        ret.fnames_r2[i] = strdup(tmp_buffer);
        ret.tmp_out_handles_r1[i] = fopen(ret.fnames_r1[i], "w");
        ret.tmp_out_handles_r2[i] = fopen(ret.fnames_r2[i], "w");
    }
    return ret;
}


static inline void splitmark_core(kseq_t *seq1, kseq_t *seq2, kseq_t *seq_index,
				    mss_settings_t settings, mark_splitter_t splitter)
{
	int l1, l2, l_index, bin;
	int count = 0;
	int pass;
    while ((l1 = kseq_read(seq1)) >= 0) {
        count += 1;
        if(!(count % settings.notification_interval)) {
            fprintf(stderr, "Number of records processed: %i.\n", count);
        }
        // Iterate through second fastq file.
        l_index = kseq_read(seq_index);
        l2 = kseq_read(seq2);
        if (l_index < 0) {
            fprintf(stderr, "Index return value for kseq_read less than "
                            "0. Are the fastqs different sizes? Abort!\n");
            FREE_SETTINGS(settings);
            FREE_SPLITTER(splitter);
            exit(EXIT_FAILURE);
        }
        if (l2 < 0) {
            fprintf(stderr, "Read 2 return value for kseq_read less than "
                            "0. Are the fastqs different sizes? Abort!\n");
            FREE_SETTINGS(settings);
            FREE_SPLITTER(splitter);
            exit(EXIT_FAILURE);
        }
        pass = test_hp(seq_index, settings.hp_threshold);
        //fprintf(stdout, "Randomly testing to see if the reading is working. %s", seq1->seq.s);
        bin = get_binner(seq_index->seq.s, settings.n_nucs);
        KSEQ_2_FQ(splitter.tmp_out_handles_r1[bin], seq1, seq_index, pass);
        KSEQ_2_FQ(splitter.tmp_out_handles_r2[bin], seq2, seq_index, pass);
    }
}
