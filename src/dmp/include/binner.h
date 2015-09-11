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


#define char_to_num(character, increment) switch(character) {\
        case 'C' : increment = 1; break;\
        case 'G' : increment = 2; break;\
        case 'T' : increment = 3; break;\
        default: increment = 0; break;\
    }


inline uint64_t get_binnerl(char *barcode, int length)
{
    uint64_t bin = 0;
    int inc_binner;
    size_t count = 0;
    for(int i = length;i > 0;i--){
        char_to_num(barcode[i - 1], inc_binner);
        bin += ((count << 2) * inc_binner);
        count++;
    }
    return bin;
}


inline int get_binner(char *barcode, int length)
{
    fprintf(stderr, "Get bin.\n");
	char *omgz = (char *)malloc((length + 1) * sizeof(char));
	omgz[length] = '\0';
	memcpy(omgz, barcode, length);
    int bin = 0;
    int inc_binner;
    size_t count = 1;
    for(int i = length; i; i--){
        char_to_num(barcode[i - 1], inc_binner);
        fprintf(stderr, "inc_binner is %i.\n", inc_binner);
        bin += ( (count << 2) * inc_binner);
        count++;
    }
    fprintf(stderr, "substr: %s. length: %i.\n", omgz, length);
    free(omgz);
    return bin;
}
