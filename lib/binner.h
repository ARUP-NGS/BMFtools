#ifndef BINNER_H
#define BINNER_H
#include <stdint.h>
#include <inttypes.h>
#include "dlib/char_util.h"
#include "dlib/compiler_util.h"
#include "dlib/math_util.h"

/* get_binner is written in a type-generic way.
 * You must declare the binner with DECLARE_BINNER and then use
 * get_binner_type to access the correct function.
 */
#define get_binner_type(barcode, length, type_t) get_binner_##type_t(barcode, length)
// get_binner defaults to uint64_t for its type.
#define get_binner(barcode, length) get_binner_uint64_t(barcode, length)

#define DECLARE_BINNER(type_t) \
    CONST static inline type_t get_binner_##type_t(char *barcode, size_t length) {\
        type_t bin = 0;\
        barcode += length;\
		LOG_DEBUG("ipow ret: %i\n", (int)ipow(4, length));
        while(length--) bin += ipow(4, length) * nuc2num_acgt(*--barcode);\
		LOG_DEBUG("bin: %lu.\n", bin);\
        return bin;\
    }

DECLARE_BINNER(uint64_t)

#endif /* BINNER_H */
