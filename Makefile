######################################
# Makefile written by Zev Kronenberg #
#     zev.kronenberg@gmail.com       #
######################################

CC=gcc
GIT_VERSION := $(shell git describe --abbrev=4 --dirty --always)
FLAGS=-std=gnu11 -Wall -fopenmp -DVERSION=\"$(GIT_VERSION)\" #-fpermissive # -fpermissive is crazy
LD= -lm -lz -lpthread 
INCLUDE= -Isrc -Ihtslib -Ihtslib/htslib -I. -Ilib -Iinclude -Iinclude/cephes
LIB=-Lhtslib

OPT_FLAGS = -O3 -DNDEBUG -finline-functions
DB_FLAGS = -fno-inline
GP_FLAGS = $(DB_FLAGS) -pg

.PHONY: all clean

all: lh3sort libhts.a fqmarksplit crms bmfsort bin

lh3sort:
	cd include/sort && make && cd ../..

libhts.a:
	cd htslib && make && cp libhts.a ../
fqmarksplit:
	$(CC) $(FLAGS) $(INCLUDE) $(LIB) $(LD) $(DB_FLAGS) src/fqmarksplit.c src/uthash_dmp_core.c -o fqmarksplit
	$(CC) $(FLAGS) $(INCLUDE) $(LIB) $(LD) $(DB_FLAGS) src/fqmarksplit.c src/uthash_dmp_core.c -o fqmarksplit_db
	$(CC) $(FLAGS) $(INCLUDE) $(LIB) $(LD) $(GP_FLAGS) src/fqmarksplit.c src/uthash_dmp_core.c -o fqmarksplit_p
crms:
	$(CC) $(FLAGS) $(INCLUDE) $(LIB) $(LD) $(OPT_FLAGS) src/crms.c src/igamc_cephes.c uthash_dmp_core.c -o crms
	$(CC) $(FLAGS) $(INCLUDE) $(LIB) $(LD) $(DB_FLAGS) src/crms.c src/igamc_cephes.c uthash_dmp_core.c -o crms_db
	$(CC) $(FLAGS) $(INCLUDE) $(LIB) $(LD) $(GP_FLAGS) src/crms.c src/igamc_cephes.c uthash_dmp_core.c -o crms_p
bmfsort:
	$(CC) $(FLAGS) $(INCLUDE) $(LIB) $(LD) $(OPT_FLAGS) src/bmfsort.c -o bmfsort
	$(CC) $(FLAGS) $(INCLUDE) $(LIB) $(LD) $(DB_FLAGS) src/bmfsort.c -o bmfsort_db
	$(CC) $(FLAGS) $(INCLUDE) $(LIB) $(LD) $(GP_FLAGS) src/bmfsort.c -o bmfsort_p
bin:
	mv crms crms_db crms_p bmfsort bmfsort_db bmfsort_p lh3sort fqmarksplit fqmarksplit_db fqmarksplit_p bin/


clean:
	rm *~ && rm *.a && rm *.o && cd htslib && make clean
