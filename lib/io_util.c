#include "io_util.h"

int isfile(char *fname)
{
	return access(fname, F_OK) != -1;

}
