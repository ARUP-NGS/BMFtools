#ifndef CSTR_UTIL_H
#define CSTR_UTIL_H
#include <inttypes.h>
#include <time.h>
#include <zlib.h>
#include "kseq.h"
#include "stdio.h"
#include "stddef.h"
#include "stdint.h"
#include "math.h"
#include "charcmp.h"

#ifndef KSEQ_DEC_GZ
#define KSEQ_DEC_GZ
KSEQ_INIT(gzFile, gzread)
#endif

#ifndef MAX_BARCODE_LENGTH
#define MAX_BARCODE_LENGTH 30
#endif

#ifndef SEQBUF_SIZE
#define SEQBUF_SIZE 300
#endif

static char *rand_string(char *str, size_t size)
{
	srand(time(NULL)); // Pick a seed!
    const char charset[] = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKSTFUOMGZWTF";
    if (size) {
        --size;
        for (size_t n = 0; n < size; n++) {
            str[n] = charset[(int)(rand() % (int) (sizeof charset - 1))];
        }
        str[size] = '\0';
    }
    return str;
}

static void fill_csv_buffer(int readlen, uint32_t *arr, char *buffer, char *prefix_typecode)
{
	char tmpbuf[12];
	strcpy(buffer, prefix_typecode);
	for(uint32_t i = 0; i < readlen; i++) {
		sprintf(tmpbuf, ",%i", arr[i]);
		strcat(buffer, tmpbuf);
	}
}


static inline void append_csv_buffer(int readlen, uint32_t *arr, char *buffer, char *prefix_typecode)
{
	char tmpbuf[12];
	strcat(buffer, prefix_typecode);
	for(int i = 0; i < readlen; i++) {
		sprintf(tmpbuf, ",%" PRIu32 "", arr[i]);
		strcat(buffer, tmpbuf);
	}
}

static inline void append_int_tag(char *buffer, char *tag, int i)
{
    char tmpbuf[15];
	sprintf(tmpbuf, "\t%s:i:%i", tag, i);
    strcat(buffer, tmpbuf);
}


static inline void fill_pv(int readlen, uint32_t *arr, char *buffer)
{
	return fill_csv_buffer(readlen, arr, buffer, (char *)"PV:B:I");
}

static inline void fill_fa(int readlen, uint32_t *arr, char *buffer)
{
	return fill_csv_buffer(readlen, arr, buffer, (char *)"FA:B:I");
}


static inline void kfill_rc(kseq_t *seq, char *buffer) {
	uint32_t i;
	for(i = 0; i < seq->seq.l; ++i){
		buffer[i] = nuc_cmpl(seq->seq.s[seq->seq.l - i - 1]);
	}
	memcpy(buffer + seq->seq.l, (char *)"\n+\n", 3 * sizeof(char));
	for(i = 0; i < seq->qual.l; ++i) {
		buffer[i + seq->qual.l + 3] = seq->qual.s[seq->qual.l - i - 1];
	}
	buffer[(seq->seq.l << 1) + 3] = '\0';
}


static inline void fill_rc(char *str, char *buffer, int len) {
#if !NDEBUG
	fprintf(stderr, "Now filling buffer at %p with string's revcmp %s.\n", &buffer[0], str);
#endif
	for(uint32_t i = 0; str[i]; ++i){
		buffer[i] = nuc_cmpl(str[len - i - 1]);
	}
	buffer[len] = '\0';
}

static inline void fill_rv(char *str, char *buffer, int len) {
	for(uint32_t i = 0; str[i]; ++i){
		buffer[i] = str[len - i - 1];
	}
	buffer[len] = '\0';
}

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


static inline int lex_memcmp(char *s1, char *s2, size_t l)
{
	for(uint64_t i = 0; i < l; ++i) {
		if(s1[i] < s2[i]) {
			return 1;
		}
		else if(s2[i] < s1[i]) {
			return 0;
		}
	}
	return -1;
}

static inline int lex_strlt(char *s1, char *s2)
{
	for(uint64_t i = 0; s1[i]; ++i) {
		if(s1[i] < s2[i]) {
			//fprintf(stderr,"String '%s' is < String '%s'. Is that right? How many fucks do I give?\n", s1, s2);
			return 0;
		}
		else if(s2[i] < s1[i]) {
			//fprintf(stderr,"String '%s' is > String '%s'. Is that right? How many fucks do I give?\n", s1, s2);
			return 1;
		}
	}
	return -1;
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

#endif
