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
