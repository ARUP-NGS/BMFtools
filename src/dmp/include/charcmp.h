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

inline void nuc_cmpl(char character, char ret) {
    switch(character) {
        case 'A': ret = 'T'; return;
        case 'C': ret = 'G'; return;
        case 'G': ret = 'C'; return;
        case 'T': ret = 'A'; return;
        default: ret = 'N'; return;
    }
}


inline int nuc_cmp(char forward, char reverse)
{
    return forward - reverse;
}
