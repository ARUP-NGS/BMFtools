#include <cstring>
#include <cstdio>
#include <stddef.h>
#include <cstdint>
#include <cstdlib>
#include "dlib/logging_util.h"

extern "C" {
    extern int sort_main(int argc, char *argv[]);
}

namespace bmf {
    extern int cap_main(int argc, char *argv[]);
    extern int depth_main(int argc, char *argv[]);
    extern int dmp_main(int argc, char *argv[]);
    extern int err_main(int argc, char *argv[]);
    extern int famstats_main(int argc, char *argv[]);
    extern int filter_main(int argc, char *argv[]);
    extern int hashdmp_main(int argc, char *argv[]);
    extern int hashdmp_inmem_main(int argc, char *argv[]);
    extern int mark_main(int argc, char *argv[]);
    extern int rsq_main(int argc, char *argv[]);
    extern int sdmp_main(int argc, char *argv[]);
    extern int stack_main(int argc, char *argv[]);
    extern int target_main(int argc, char *argv[]);
    extern int vet_main(int argc, char *argv[]);
}
