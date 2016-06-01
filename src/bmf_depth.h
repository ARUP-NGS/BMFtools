#ifndef BMF_DEPTH_H
#define BMF_DEPTH_H
#include <cstdint>
extern "C" {
    #include "htslib/khash.h"
}

#define DEFAULT_MAX_DEPTH (1 << 18)
namespace bmf {
    KHASH_MAP_INIT_INT(depth, uint64_t)
    int depth_main(int argc, char *argv[]);
}

#endif /* BMF_DEPTH_H */
