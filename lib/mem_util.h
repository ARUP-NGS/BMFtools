#ifndef MEM_UTIL_H
#define MEM_UTIL_H

#include <stdio.h>
#include <stdint.h>
#include <stddef.h>

// Leave semicolon out so that it looks like a normal function.
#if !NDEBUG
#define cond_free(var) do { if(var) {fprintf(stderr, "About to free variable at %p (%s).\n", var, #var); free(var); var = NULL;}} while(0)
#else
#define cond_free(var) do {if(var) {free(var); var = NULL;}} while(0)
#endif

#define cond_realloc_t(var, newsize, type_t)\
	do {\
		var = (type_t)realloc(var, newsize);\
		if(!var){\
			fprintf(stderr, "Could not allocate new memory for size %" PRIu64 ". Abort mission!\n", newsize);\
			exit(EXIT_FAILURE);\
		}\
	} while(0)

#define cond_realloc(var, newsize)\
	do {\
		var = realloc(var, newsize);\
		if(!var){\
			fprintf(stderr, "Could not allocate new memory for size %i. Abort mission!\n");\
			exit(EXIT_FAILURE);\
		}\
	} while(0)

#define cond_close(fp) \
	if(fp) {\
		fclose(fp);\
	}\
	fp = NULL

#define roundup_div(top, bottom) 1 + (((top) - 1) / (bottom))


#endif
