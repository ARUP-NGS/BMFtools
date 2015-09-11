#include "stdint.h"


// Functions
inline int ipow(int base, int exp)
{
    int result = 1;
    while (exp)
    {
        if (exp & 1)
            result *= base;
        exp >>= 1;
        base *= base;
    }

    return result;
}


inline void char_to_num(char character, int increment) {
    switch(character) {
        case 'C' : increment = 1; return;
        case 'G' : increment = 2; return;
        case 'T' : increment = 3; return;
        default: increment = 0; return;
    }
}


inline uint64_t get_binnerl(char *barcode, int length)
{
    uint64_t bin = 0;
    int inc_binner;
    size_t count = 0;
    for(int i = length;i > 0;i--){
        char_to_num(barcode[i - 1], inc_binner);
        bin += (ipow(4, count) * inc_binner);
        count++;
    }
    return bin;
}


inline int get_binner(char *barcode, int length)
{
    int bin = 0;
    int inc_binner;
    size_t count = 0;
    for(int i = length;i > 0;i--){
        char_to_num(barcode[i - 1], inc_binner);
        bin += (ipow(4, count) * inc_binner);
        count++;
    }
    return bin;
}
