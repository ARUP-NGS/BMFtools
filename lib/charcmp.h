#ifndef CHARCMP_H
#define CHARCMP_H

static const uint32_t nucpos_arr[128] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
										 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,2,0,0,0,0,0,0,
										 0,0,0,0,0,0,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,3,
										 0,0,0,0,0,0,0,0,0,0,0};

static const char *rc_string = "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTNGNNNCNNNNNNNNNNNNANNNNNNNNNNNNTNGNNNCNNNNNNNNNNNNANNNNNNNNNNN";

#define nuc2num(character) nucpos_arr[(int)character]

#define nuc_cmpl(character) rc_string[(int)character]


static inline int nuc_cmp(char forward, char reverse)
{
	return forward - reverse;
}

static inline uint32_t nuc_to_posdata(char character)
{
	return nuc2num(character) | ((character & 0x8U) >> 1); // This sets bit 3 (4) to true if the base is N or n, out of ACGTNacgtn
}


//Converts a nucleotide in a char * into an index for the phred_sums and nuc_counts arrays.
static inline void nuc_to_pos(char character, int *nuc_indices)
{
	switch(character) {
		case 'A': nuc_indices[0] = 0; nuc_indices[1] = 0; return;
		case 'C': nuc_indices[0] = 1; nuc_indices[1] = 1; return;
		case 'G': nuc_indices[0] = 2; nuc_indices[1] = 2; return;
		case 'T': nuc_indices[0] = 3; nuc_indices[1] = 3; return;
		default: nuc_indices[0] = 0; nuc_indices[1] = 4; return;
	}
}

#endif
