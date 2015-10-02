#include <stdio.h>
#include <stdint.h>
#include <stddef.h>

// Leave semicolon out so that it looks like a normal function.
#ifndef cond_free
#define cond_free(var) do {if(var) {free(var); var = NULL;}} while(0)
#endif
