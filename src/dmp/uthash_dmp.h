#include "include/kseq_dec.h"

/*
 * :param: seq - [arg/kseq_t *] a filled-in kseq object.
 * :param: buf - a pre-allocated buffer or malloc'd char_ptr with enough space for the barcode and the null terminus.
 * :returns:
 */
inline void cp_bs2buf(kseq_t *seq, char *buf)
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
