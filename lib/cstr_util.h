#ifndef CSTR_UTIL_H
#define CSTR_UTIL_H
#include "kseq.h"

#ifndef KSEQ_DEC_GZ
#define KSEQ_DEC_GZ
KSEQ_INIT(gzFile, gzread)
#endif


static inline void fill_csv_buffer(int readlen, int *arr, char *buffer, char *prefix_typecode)
{
	char tmpbuf[20];
	strcpy(buffer, prefix_typecode);
	for(int i = 0; i < readlen; i++) {
		sprintf(tmpbuf, ",%i", arr[i]);
		strcat(buffer, tmpbuf);
	}
}

static inline char nuc_cmpl(char character) {
	switch(character) {
		case 'A': return 'T';
		case 'C': return 'G';
		case 'G': return 'C';
		case 'T': return 'A';
		default: return character;
	}
}

static inline void append_csv_buffer(int readlen, uint32_t *arr, char *buffer, char *prefix_typecode)
{
	char tmpbuf[20];
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
	return fill_csv_buffer(readlen, (int *)arr, buffer, (char *)"PV:B:I");
}

static inline void fill_fa(int readlen, uint32_t *arr, char *buffer)
{
	return fill_csv_buffer(readlen, (int *)arr, buffer, (char *)"FA:B:I");
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

#endif
