#pragma once
#include "stdio.h"
#include "stddef.h"
#include "stdint.h"
#include "math.h"
#include "charcmp.h"
#include "khash.h"
#include "uthash.h"

/*
 * Returns a null-terminated string with the extension and terminal period removed.
 * Warning: Must be freed!
 */
inline char *trim_ext(char *fname) {
    char *buf = malloc((strlen(fname) + 1) * sizeof(char ));
    ptrdiff_t pos = strrchr(fname, '.') - fname; // Find the position in the read where the last '.' is.
    memcpy(buf, fname, pos * sizeof(char));
    buf[pos] = '\0';
    fprintf(stderr, "tmp buffer: %s.\n", buf);
    return buf;
}
