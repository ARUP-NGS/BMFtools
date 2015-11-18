#include "crms.h"

static void print_usage(char *argv[]);
static void print_opt_err(char *argv[], char optarg[]);
static double igamc_pvalues(int num_pvalues, double x);
KingFisher_t init_kf(int readlen);
void pushback_kseq(KingFisher_t *fisher, kseq_t *seq, int *nuc_indices, int blen);
static int bmftools_dmp_core(kseq_t *seq, FILE *out_handle);
int ARRG_MAX(KingFisher_t *kfp, int index);
char ARRG_MAX_TO_NUC(int argmaxret);
int infer_barcode_length(char *bs_ptr);
void destroy_kf(KingFisher_t *kfp);
void clear_kf(KingFisher_t *kfp);
int get_binner(char *barcode, int length);
void nuc_to_pos(char character, int *nuc_indices);
char test_hp(char *seq, int threshold);
char rescale_qscore(int readnum, int qscore, int cycle, char base, int readlen, char *rescaler);

static int bmftools_dmp_wrapper(char *input_path, char *output_path);
