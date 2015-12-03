#ifndef KHASH_DMP_CORE_H
#define KHASH_DMP_CORE_H
#include "khash.h"
#include "crms.h"

KHASH_MAP_INIT_STR(dmp, KingFisher_t *)
void khash_dmp_core(char *infname, char *outfname);
int khash_dmp_main(int argc, char *argv[]);
void splitterhash_destroy(splitterhash_params_t *params);
splitterhash_params_t *init_splitterhash(mssi_settings_t *settings_ptr, mark_splitter_t *splitter_ptr);
splitterhash_params_t *init_splitterhash_mss(mss_settings_t *settings_ptr, mark_splitter_t *splitter_ptr);

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


static inline tmpvars_t *init_tmpvars_p(char *bs_ptr, int blen, int readlen)
{
	tmpvars_t *ret = (tmpvars_t *)malloc(sizeof(tmpvars_t));
	ret->blen = blen;
	ret->readlen = readlen;
	ret->bs_ptr = bs_ptr;
	ret->buffers = (tmpbuffers_t *)malloc(sizeof(tmpbuffers_t));
	ret->buffers->name_buffer[0] = '@';
	ret->buffers->name_buffer[blen] = '\0';
	return ret;
}


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
	char *view = barcode_mem_view(seq);
	int blen = 0;
	while(view[blen] != '\0' && view[blen] != '|') {
		++blen;
	}
	memcpy(buf, view, blen);
	buf[blen] = '\0';
	return;
}



#endif // KHASH_DMP_CORE_H
