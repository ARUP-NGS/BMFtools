#ifndef BINNER_H
#define BINNER_H
#include "stdint.h"
#include "inttypes.h"


static inline uint64_t ulpow(uint64_t base, uint64_t exp);
static inline int ipow(int base, int exp);
// Functions
static inline int ipow(int base, int exp)
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

/* get_binner is written in a type-generic way.
 * You must declare the binner with DECLARE_BINNER and then use
 * get_binner_type to access the correct function.
 */
#define get_binner_type(barcode, length, type_t) get_binner_##type_t(barcode, length)

#define DECLARE_BINNER(type_t) \
	CONST static inline type_t get_binner_##type_t(char *barcode, size_t length) {\
		type_t bin = 0;\
		uint64_t i;\
		for(i = 0; i < length; ++i)\
			bin += ulpow(4, i) * nuc2num_acgt(*barcode++);\
		return bin;\
	}

DECLARE_BINNER(uint64_t)


static inline int64_t lpow(int64_t base, int64_t exp)
{
	int64_t result = 1;
	while (exp)
	{
		if (exp & 1)
			result *= base;
		exp >>= 1;
		base *= base;
	}

	return result;
}


static inline uint64_t ulpow(uint64_t base, uint64_t exp)
{
	uint64_t result = 1;
	while (exp)
	{
		if (exp & 1)
			result *= base;
		exp >>= 1;
		base *= base;
	}
	return result;
}

#endif /* BINNER_H */
