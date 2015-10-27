#ifndef CSTR_UTILS_H
#define CSTR_UTILS_H

#include <inttypes.h>
#include "stdio.h"
#include "stddef.h"
#include "stdint.h"
#include "math.h"
#include "charcmp.h"

#ifndef MAX_BARCODE_LENGTH
#define MAX_BARCODE_LENGTH 30
#endif

/*
 * Returns a null-terminated string with the extension and terminal period removed.
 * Warning: Must be freed!
 */
static inline char *trim_ext(char *fname)
{
	fprintf(stderr, "Now trimming char * %s.\n", fname);
	char *buf = malloc((strlen(fname) + 1) * sizeof(char ));
	char *found_pos = strrchr(fname, '.');
	if(!found_pos) {
		fprintf(stderr, "Could not trim file name's extension. Looks like it's missing a '.' (name: '%s').\n", fname);
		exit(EXIT_FAILURE);
	}
	ptrdiff_t pos = found_pos - fname; // Find the position in the read where the last '.' is.
	memcpy(buf, fname, pos * sizeof(char));
	buf[pos] = '\0';
	//fprintf(stderr, "tmp buffer: %s.\n", buf);
	return buf;
}

static inline void rc_barcode(char *barcode, int blen) {
	char buf[MAX_BARCODE_LENGTH];
	memcpy(buf, barcode, blen * sizeof(char));
	buf[blen] = '\0';
	for(int i = 0; i < blen; ++i) {
		barcode[i] = nuc_cmpl(buf[blen - i - 1]);
	}
	return;
}
/*
 * Fast positive atoi
 */
static inline int fp_atoi(char *str)
{
	int ret = *str++ - '0';
	while(*str) {
		ret = ret*10 + (*str++ - '0');
	}
	return ret;
}

static const char digit_pairs[201] = {
  "00010203040506070809"
  "10111213141516171819"
  "20212223242526272829"
  "30313233343536373839"
  "40414243444546474849"
  "50515253545556575859"
  "60616263646566676869"
  "70717273747576777879"
  "80818283848586878889"
  "90919293949596979899"
};


static inline char *opt_itoa(unsigned val, char *s)
{
	if(!val)
	{
		s="0";
		return s;
	}

	int size;
	if(val>=10000)
	{
		if(val>=10000000)
		{
			if(val>=1000000000)
				size=10;
			else if(val>=100000000)
				size=9;
			else
				size=8;
		}
		else
		{
			if(val>=1000000)
				size=7;
			else if(val>=100000)
				size=6;
			else
				size=5;
		}
	}
	else
	{
		if(val>=100)
		{
			if(val>=1000)
				size=4;
			else
				size=3;
		}
		else
		{
			if(val>=10)
				size=2;
			else
				size=1;
		}
	}

	s[size] = '\0';
	char* c = &s[size-1];
	while(val>=100)
	{
	   int pos = val % 100;
	   val /= 100;
	   *(short*)(c-1)=*(short*)(digit_pairs+2*pos);
	   c-=2;
	}
	while(val>0)
	{
		*c--='0' + (val % 10);
		val /= 10;
	}
	return s;
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

static inline char *revcmp(char *dest, char *src, uint64_t l)
{
	//char *ret = malloc((l + 1)* sizeof(char));
	dest[l] = '\0';
	for(uint64_t i = 0; i < l; ++i) {
		dest[i] = nuc_cmpl(src[l - i - 1]);
	}
	return dest;
}

static inline int lex_lt(char *s, size_t l)
{
	//fprintf(stderr, "Char *: %s.\n", s);
	for(uint64_t i = 0; i < l; ++i) {
		//fprintf(stderr, "Comparing string indices %"PRIu64" and %"PRIu64".\n", i, l - i - 1);
		if(s[l - i - 1] > s[i]) {
			//fprintf(stderr, "Character %c (%i) at index %i is less than %c. String: %s.\n", s[i], s[i], i, s[l -i - 1], s);
			return 1;
		}
		else if(s[l - i - 1] < s[i]) {
			//fprintf(stderr, "Character %c (%i) is less than %c. String: %s.\n", s[l - i - 1], s[l - i - 1], s[i], s);
			return 0;
		}
	}
	//fprintf(stderr, "This barcode is palindromic!\n");
	return -1; // Palindromic
}

static inline int rclex_lt(char *s, size_t l)
{
	char cmp;
	//fprintf(stderr, "Char *: %s.\n", s);
	for(uint64_t i = 0; i < l; ++i) {
		cmp = nuc_cmpl(l - i - 1);
		//fprintf(stderr, "Comparing string indices %"PRIu64" and %"PRIu64".\n", i, l - i - 1);
		if(cmp > s[i]) {
			//fprintf(stderr, "Character %c (%i) at index %i is less than %c. String: %s.\n", s[i], s[i], i, s[l -i - 1], s);
			return 1;
		}
		else if(s[i] > cmp) {
			//fprintf(stderr, "Character %c (%i) is less than %c. String: %s.\n", s[l - i - 1], s[l - i - 1], s[i], s);
			return 0;
		}
	}
	//fprintf(stderr, "This barcode is palindromic!\n");
	return -1; // RC palindromic
}

#endif
