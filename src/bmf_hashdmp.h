#ifndef BMF_HASHDMP_H
#define BMF_HASHDMP_H
#include <assert.h>

#include "lib/mseq.h"
#ifdef __cplusplus
extern "C" {
#endif
#include "include/uthash.h"
#include "lib/kingfisher.h"
#ifdef __cplusplus
}
#endif

KHASH_MAP_INIT_STR(dmp, KingFisher_t *)
void hash_dmp_core(char *infname, char *outfname);
int hash_dmp_main(int argc, char *argv[]);
#ifdef __cplusplus
extern "C" {
#endif
	extern void splitterhash_destroy(splitterhash_params_t *params);
	extern splitterhash_params_t *init_splitterhash(marksplit_settings_t *settings_ptr, mark_splitter_t *splitter_ptr);
	void stranded_hash_dmp_core(char *infname, char *outfname);
	tmpvars_t *init_tmpvars_p(char *bs_ptr, int blen, int readlen);
#ifdef __cplusplus
}
#endif

typedef struct HashKing {
	UT_hash_handle hh;
	char id[MAX_BARCODE_LENGTH + 1];
	KingFisher_t *value;
}hk_t;


CONST static inline int infer_barcode_length(char *bs_ptr)
{
	char *const current = bs_ptr;
	for (;;) {
		switch(*bs_ptr++) {
		case '|': // Fall-through
		case '\0': return bs_ptr - current - 1;
		}
	}
	return -1; // This never happens.
}

static inline void cp_view2buf(char *view, char *buf)
{
	for(;;) {
		switch(*view) {
			case '\0': // Fall-through
			case '|': *buf++ = '\0'; return;
			default: *buf++ = *view++;
		}
	}
}

static inline KingFisher_t *init_kfp(size_t readlen)
{
	const size_t r5 = readlen * 5;
	KingFisher_t *ret = (KingFisher_t *)calloc(1, sizeof(KingFisher_t));
	ret->readlen = readlen;
	ret->max_phreds = (char *)malloc((r5) * sizeof(char));
	ret->nuc_counts = (uint16_t *)calloc(r5, sizeof(uint16_t));
	ret->phred_sums = (uint32_t *)calloc(r5, sizeof(uint32_t));
	memset(ret->max_phreds, '#', r5);
	ret->pass_fail = '1';
	return ret;
}

static inline void destroy_kf(KingFisher_t *kfp)
{
	free(kfp->nuc_counts);
	free(kfp->phred_sums);
	free(kfp->max_phreds);
	free(kfp);
}

tmpvars_t *init_tmpvars_p(char *bs_ptr, int blen, int readlen);


static inline void tmpvars_destroy(tmpvars_t *tmp)
{
	free(tmp->buffers), free(tmp);
}


/* @func
 * Copies the barcode sequence from a fastq comment field into a buffer
 * :param: seq - [arg/kseq_t *] a filled-in kseq object.
 * :param: buf - a pre-allocated buffer or malloc'd char_ptr with enough space for the barcode and the null terminus.
 */
#define cp_bs2buf(seq, buf) cp_view2buf(barcode_mem_view(seq), buf)



#endif /* BMF_HASHDMP_H */
