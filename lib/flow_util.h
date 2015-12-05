#ifndef FLOW_UTIL_H
#define FLOW_UTIL_H

#include "stdlib.h"
#include "stdio.h"

void log_abort(char *str) {
	fprintf(stderr, str);
	exit(EXIT_FAILURE);
}
#endif
