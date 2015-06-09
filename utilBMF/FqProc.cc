/*
FqMerge.cxx

simple C function that merges a fastq family.

    used to see how we can do this with Cython/numpy

*/
#include <stdlib.h>
#include "FqProc.h"


// Argmax - finds the index for which arr is maximized.
int amax(int* arr, int length){
	int max = 0;
	int amax = 0;
	for(int i = 0; i < length; i++){
		if(arr[i] > max){max = arr[i]; amax = i;}
	}
	return amax;
}

Read FqMerge (int** sequence, int** quality, int nreads, int lenreads) {
	int* newSeq = (int*)malloc(lenreads*sizeof(int));
	int* newQual = (int*)malloc(lenreads*sizeof(int));
	int* tmpQual = (int*)malloc(4*sizeof(int));  // Holds the summed phred qualities for differing nucleotide calls.
	for(int i = 0; i < lenreads; i++){
		for(int j = 0; j < nreads; j++){
			for(int k = 0; k < 4; k++){
				switch (sequence[i][j]) {
				case 65:
					tmpQual[0] += quality[i][j];
					break;
				case 67:
					tmpQual[1] += quality[i][j];
					break;
				case 71:
					tmpQual[2] += quality[i][j];
					break;
				case 84:
					tmpQual[3] += quality[i][j];
					break;
				}
			}
		}
		switch(amax(tmpQual, 4)){
		case 0:
			newSeq[i] = 65;
			newQual[i] = tmpQual[0];
			break;
		case 1:
			newSeq[i] = 67;
			newQual[i] = tmpQual[1];
			break;
		case 2:
			newSeq[i] = 71;
			newQual[i] = tmpQual[2];
			break;
		case 3:
			newSeq[i] = 84;
			newQual[i] = tmpQual[3];
			break;
		}
	}
	free(tmpQual);
	Read retRead = {newSeq, newQual, lenreads};
    return retRead;
}

int main(int argv, char** args){
	return 0;
}
