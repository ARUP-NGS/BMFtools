#include "call_lh3_sort.h"

// Set global parameters as needed.

/* Lists of key field comparisons to be tried. */
static struct keyfield keyhead;

int parallelism_enabled = 1;


void init_lh3sort_key() {
    struct keyfield *key = (struct keyfield *) xmalloc(sizeof(struct keyfield));
    key_init(key);
    key->sword = 2;
    key->sword = 3;
    key->schar = 0;
    key->echar = 0;
    key->skipsblanks = 0;
    key->skipeblanks = 0;
    key->translate = NULL;
    key->next = NULL;
    keyhead.sword = 0,
    keyhead.schar = 0,
    keyhead.eword = 0,
    keyhead.echar = 0,
    keyhead.skipsblanks = 0,
    keyhead.skipeblanks = 0,
    keyhead.translate = NULL;
    keyhead.next = key;
}

void apply_parallel_lh3sort(char **filenames, int nfilenames,
										  FILE **ofps, int nthreads);
