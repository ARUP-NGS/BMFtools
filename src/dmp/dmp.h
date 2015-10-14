#include "dmp_interface.h"

static double igamc_pvalues(int num_pvalues, double x);
KingFisher_t init_kf(int readlen);
void pushback_kseq(KingFisher_t *fisher, kseq_t *seq, int *nuc_indices, int blen);
int bmftools_dmp_core(kseq_t *seq, FILE *out_handle);
int ARRG_MAX(KingFisher_t *kfp, int index);
char ARRG_MAX_TO_NUC(int argmaxret);
int pvalue_to_phred(double pvalue);
void fill_fa_buffer(KingFisher_t *kfp, int *agrees, char *buffer);
void dmp_process_write(KingFisher_t *kfp, FILE *handle, int blen, tmpbuffers_t *tmp);
int bmftools_dmp_wrapper(char *input_path, char *output_path);
int infer_barcode_length(char *bs_ptr);
void fill_csv_buffer(int readlen, int *arr, char *buffer, char *prefix, char typecode);
void fill_fa_buffer(KingFisher_t *kfp, int *agrees, char *buffer);
void fill_pv_buffer(KingFisher_t *kfp, int *agrees, char *buffer);
void destroy_kf(KingFisher_t *kfp);
void clear_kf(KingFisher_t *kfp);
void u32toa_branchlut(uint32_t value, char* buffer);
void i32toa_branchlut(int32_t value, char* buffer);
int get_binner(char *barcode, int length);
void nuc_to_pos(char character, int *nuc_indices);
char test_hp(char *seq, int threshold);
mseq_t *p7_mseq_rescale_init(kseq_t *seq, char *rescaler, int n_len, int is_read2);
char rescale_qscore(int readnum, int qscore, int cycle, char base, int readlen, char *rescaler);
