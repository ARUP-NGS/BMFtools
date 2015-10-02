#include "uthash_dmp_core.h"
#include "include/o_mem.h"

#ifndef MAX_N_BLENS
#define MAX_N_BLENS 5
#endif
#ifndef MAX_HOMING_SEQUENCE
#define MAX_HOMING_SEQUENCE 8
#endif

typedef struct blens {
    int max_blen; // Last value in blens
    int min_blen; // Lowest value in blens
    int blens[MAX_N_BLENS]; // Array holding blens
    int n; // Number of blens to look for
    int current_blen;
    int homing_sequence_length;
    char homing_sequence[MAX_HOMING_SEQUENCE + 1];
} blens_t;


typedef struct crms_settings {
    int hp_threshold; // The minimum length of a homopolymer run to fail a barcode.
    int n_nucs; // Number of nucleotides to split by.
    char *output_basename;
    char *input_r1_path;
    char *input_r2_path;
    int n_handles; // Number of handles
    int notification_interval; // How many sets of records do you want to process between progress reports?
    blens_t *blen_data;
    int offset; // Number of bases at the start of the inline barcodes to skip for low quality.
    char ****rescaler; // Three-dimensional rescaler array. Size: [readlen, 39, 4] (length of reads, number of original quality scores, number of bases).
    char *rescaler_path; // Path to flat text file for parsing in the rescaler.
    char *ffq_prefix; // Final fastq prefix.
    int threads; // Number of threads to use for parallel dmp.
    char *homing_sequence;
    int homing_sequence_length;
    int run_hash_dmp;
} crms_settings_t;
/*
typedef struct mssi_settings {
    char *homing_sequence; // Homing sequence...
    int homing_sequence_length; // Length of homing sequence, should it be used.
    int threads;
    int run_hash_dmp;
} mssi_settings_t;
*/

/*
 * :param: settings [crms_settings_t, mssi_settings_t] Settings struct in which to free the rescaler.
 * :return: void
 * This function supersedes free_rescaler_array by being type-generic.
 */

#define cfree_rescaler(settings) \
    do {\
        if(settings.rescaler) {\
            int readlen##_settings = count_lines(settings.rescaler_path);\
            for(int i = 0; i < 2; ++i) {\
                for(int j = 0; j < readlen##_settings; ++j) {\
                    for(int k = 0; k < 39; ++k) {\
                        cond_free(settings.rescaler[i][j][k]);\
                    }\
                    cond_free(settings.rescaler[i][j]);\
                }\
                cond_free(settings.rescaler[i]);\
            }\
            cond_free(settings.rescaler);\
        }\
    } while(0)

void free_rescaler_array(crms_settings_t settings) {
    int readlen = count_lines(settings.rescaler_path);
    for(int i = 0; i < 2; ++i) {
        for(int j = 0; j < readlen; ++j) {
            for(int k = 0; k < 39; ++k) {
                cond_free(settings.rescaler[i][j][k]);
            }
            cond_free(settings.rescaler[i][j]);
        }
        cond_free(settings.rescaler[i]);
    }
    cond_free(settings.rescaler);
    return;
}

// Inline function declarations
int crc_flip(mseq_t *mvar, char *barcode, int blen, int readlen);
void crc_mseq(mseq_t *mvar, tmp_mseq_t *tmp);
void FREE_SPLITTER(mark_splitter_t var);
mseq_t init_rescale_revcmp_mseq(kseq_t *seq, char *barcode, char ****rescaler, tmp_mseq_t *tmp, int n_len, int is_read2);
mark_splitter_t init_splitter_inline(mssi_settings_t* settings_ptr);
mark_splitter_t init_splitter(mss_settings_t *settings_ptr);
int ipow(int base, int exp);
void mseq_destroy(mseq_t *mvar);
void mseq_rescale_init(kseq_t *seq, mseq_t *ret, char ****rescaler, tmp_mseq_t *tmp, int n_len, int is_read2);
int nuc2num(char character);
int nuc_cmp(char forward, char reverse);
char ****parse_rescaler(char *qual_rescale_fname);
int rescale_qscore(int readnum, int qscore, int cycle, char base, char ****rescaler);
void set_barcode(kseq_t *seq1, kseq_t *seq2, char *barcode, int offset, int blen1_2);
int test_homing_seq(kseq_t *seq1, kseq_t *seq2, mssi_settings_t *settings_ptr);
char test_hp_inline(char *barcode, int length, int threshold);
char test_hp(kseq_t *seq, int threshold);
void tmp_mseq_destroy(tmp_mseq_t mvar);
void update_mseq(mseq_t *mvar, char *barcode, kseq_t *seq, char ****rescaler, tmp_mseq_t *tmp, int n_len, int is_read2);
char nuc_cmpl(char character);
void mseq2fq_inline(FILE *handle, mseq_t *mvar, char pass_fail);
int count_lines(char *fname);
void FREE_SPLITTER(mark_splitter_t var);
char *trim_ext(char *fname);
char *make_default_outfname(char *fname, const char *suffix);
char *make_crms_outfname(char *fname);
int get_fileno_limit();
void increase_nofile_limit(int new_limit);
uint64_t get_binnerul(char *barcode, int length);
int get_binner(char *barcode, int length);
uint64_t ulpow(uint64_t base, uint64_t exp);
void splitterhash_destroy(splitterhash_params_t *params);
splitterhash_params_t *init_splitterhash(mssi_settings_t *settings_ptr, mark_splitter_t *splitter_ptr);
int vl_homing_loc(kseq_t *seq1, kseq_t *seq2, crms_settings_t *settings_ptr);
blens_t *get_blens(char *str2parse);
void free_crms_settings(crms_settings_t settings);
