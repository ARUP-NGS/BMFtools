#ifndef O_IO_UTIL_H
#define O_IO_UTIL_H
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <zlib.h>

int isfile(char *fname);

inline gzFile open_gzfile(char *infname) {
    if(strcmp(infname, "-") == 0 || strcmp(infname, "stdin") == 0) {
        fprintf(stderr, "Reading from standard in because infname is %s.\n", infname);
        return gzdopen(STDIN_FILENO, "r"); // Opens stdin.
    }
    else {
        fprintf(stderr, "Reading from %s.\n", infname);
        return gzopen(infname, "r");
    }
}

inline FILE *open_ofp(char *infname) {
    if(strcmp(infname, "-") == 0 || strcmp(infname, "stdout") == 0) {
        fprintf(stderr, "Reading from standard in because infname is %s.\n", infname);
        return stdout; // Opens stdin.
    }
    else {
        fprintf(stderr, "Reading from %s.\n", infname);
        return fopen(infname, "r");
    }
}

#endif
