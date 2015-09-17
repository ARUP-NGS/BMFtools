#pragma once
//#include "branchlut.h"
/*
 * This is broken.
inline void comma_i32toa(int32_t value, char *buffer)
{
    *buffer++ = ',';
    i32toa_branchlut(value, buffer);
}
*/
/*
 * Small character comparison or conversion utilities.
 */

inline int nuc2num(char character)
{
    switch(character) {
        case 'C': return 1; break;
        case 'G': return 2; break;
        case 'T': return 3; break;
        default: return 0; break; // 'A'
    }
}

inline char nuc_cmpl(char character) {
    switch(character) {
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
        default: return character;
    }
}


inline int nuc_cmp(char forward, char reverse)
{
    return forward - reverse;
}


//Converts a nucleotide in a char * into an index for the phred_sums and nuc_counts arrays.
inline void nuc_to_pos(char character, int *nuc_indices)
{
    switch(character) {
        case 'A': nuc_indices[0] = 0; nuc_indices[1] = 0; return;
        case 'C': nuc_indices[0] = 1; nuc_indices[1] = 1; return;
        case 'G': nuc_indices[0] = 2; nuc_indices[1] = 2; return;
        case 'T': nuc_indices[0] = 3; nuc_indices[1] = 3; return;
        default: nuc_indices[0] = 0; nuc_indices[1] = 4; return;
    }
}
