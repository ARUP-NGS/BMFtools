
// Typedefs
typedef struct fq_core {
    char * seq; // Sequence
    int32_t * qual; // Array of new quality scores
    int16_t rlen; //Read length
} fq_core_t;

typedef struct dmp_stack {
    char * comment;
    fq_core * recs;
    int n;
    int8_t ** max_qscores; // Maximum q score seen at position with given base calls.
} dmp_stack_t;

typedef struct fq_record {
    fq_core_t core;
    char * rname;
    char * comment;
    char * qualstr; // Not equivalent to qual in core!
} fq_record_t;

// Macros

#define STRING_BUFSIZE 1024

#define FQ_RECORD_TO_STRING(recptr, handle) fprintf(handle, "%s %s\n%s\n+\n%s", rec->rname, rec->comment,\
                                                    rec->core.seq, rec->qualstr);

// TODO:
// Function or macro to turn a read line into a fq_core
// Function or macro to convert a fq_record into a string to output
//
#include <zlib.h>
#include <stdio.h>
#include <stdint.h>
#include <regex.h>
#include "kseq.h"
#include "cephes_ll/igamil.c"


