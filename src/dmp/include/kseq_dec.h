#include "kseq.h"
#include <zlib.h>
/*
 * Utilities regarding kseq types.
 */
KSEQ_INIT(gzFile, gzread)



/*
 * Warning: returns a NULL upon not finding a second pipe symbol.
 * This is *NOT* a properly null-terminated string.
 */
char *barcode_mem_view(kseq_t *seq)
{
    int hits = 0;
    for(int i = 0; i < seq->comment.l; i++) {
        if(seq->comment.s[i] == '|') {
            if(!hits) {
                hits += 1;
            }
            else {
                return (char *)(seq->comment.s + i + 4); // 4 for "|BS="
            }
        }
    }
    return NULL;
}
