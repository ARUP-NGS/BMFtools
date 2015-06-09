/*
FqMerge.hxx

simple C header file that merges a fastq family.

    used to see how we can do this with Cython/numpy

*/

#include <stdlib.h>

typedef struct Read_ {
	int* sequence;
	int* quality;
	int length;
} Read;

int amax(int* arr, int length);
Read FqMerge (int** sequence, int** quality, int nreads, int lenreads);
