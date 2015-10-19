#ifndef CSTR_UTILS_H
#define CSTR_UTILS_H

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
static inline char *trim_ext(char *fname)
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

#endif
