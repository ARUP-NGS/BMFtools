#include "uthash_dmp_core.c"
#include "khash_dmp_core.h"

int uthash_dmp_main(int argc, char *argv[]);

int main(int argc, char *argv[])
{
	if(argc > 1 && strcmp(argv[1], "khash") == 0) {
		return khash_dmp_main(argc - 1, argv + 1);
	}
	return uthash_dmp_main(argc, argv);
}
