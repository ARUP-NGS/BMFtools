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
inline char *trim_ext(char *fname)
{
	fprintf(stderr, "Now trimming char * %s.\n", fname);
	char *buf = malloc((strlen(fname) + 1) * sizeof(char ));
	ptrdiff_t pos = strrchr(fname, '.') - fname; // Find the position in the read where the last '.' is.
	memcpy(buf, fname, pos * sizeof(char));
	buf[pos] = '\0';
	fprintf(stderr, "tmp buffer: %s.\n", buf);
	return buf;
}
/*
 * Fast positive atoi
 */
inline int fp_atoi(char *str)
{
	int ret = *str++ - '0';
	while(*str) {
		ret = ret*10 + (*str++ - '0');
	}
	return ret;
}


inline int fast_atoi(char *str)
{
	int ret = 0;
	int sign = 1;
	switch(*str) {
		case '-': sign = -1; break;
		case '+': break;
		default: ret = *str - '0';
	}
	++str;
	while(*str) {
		ret = ret*10 + (*str++ - '0');
	}
	return ret * sign;
}
