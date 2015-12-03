#ifndef KHASH_DMP_CORE_H
#define KHASH_DMP_CORE_H
#include "khash.h"
#include "uthash.h"
#include "crms.h"
#include <assert.h>

KHASH_MAP_INIT_STR(dmp, KingFisher_t *)
void khash_dmp_core(char *infname, char *outfname);
int khash_dmp_main(int argc, char *argv[]);
void splitterhash_destroy(splitterhash_params_t *params);
splitterhash_params_t *init_splitterhash(mssi_settings_t *settings_ptr, mark_splitter_t *splitter_ptr);
splitterhash_params_t *init_splitterhash_mss(mss_settings_t *settings_ptr, mark_splitter_t *splitter_ptr);

typedef struct HashKing {
	UT_hash_handle hh;
	char id[MAX_BARCODE_LENGTH + 1];
	KingFisher_t *value;
}hk_t;



static inline void cp_view2buf(char *view, char *buf)
{
	int blen = 0;
	while(view[blen] != '\0' && view[blen] != '|') {
		buf[blen] = view[blen];
        ++blen;
	}
	buf[blen] = '\0';
	return;
}


tmpvars_t *init_tmpvars_p(char *bs_ptr, int blen, int readlen);


static inline void tmpvars_destroy(tmpvars_t *tmp)
{
	free(tmp->buffers), free(tmp), tmp = NULL;
}


/*
 * :param: seq - [arg/kseq_t *] a filled-in kseq object.
 * :param: buf - a pre-allocated buffer or malloc'd char_ptr with enough space for the barcode and the null terminus.
 * :returns:
 */
static inline void cp_bs2buf(kseq_t *seq, char *buf)
{
	const char *view = barcode_mem_view(seq);
	int blen = 0;
	while(view[blen] != '\0' && view[blen] != '|') {
		++blen;
	}
	memcpy(buf, view, blen);
	buf[blen] = '\0';
	return;
}



#endif // KHASH_DMP_CORE_H
