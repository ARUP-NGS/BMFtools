#ifndef KHASH_DMP_CORE_H
#define KHASH_DMP_CORE_H
#include "khash.h"
#include "uthash_dmp_core.h"

KHASH_MAP_INIT_STR(dmp, KingFisher_t *)
void khash_dmp_core(char *infname, char *outfname);
int khash_dmp_main(int argc, char *argv[]);

#endif // KHASH_DMP_CORE_H
