######################################
# Makefile written by Daniel Baker   #
#     d.nephi.baker@gmail.co m       #
#  Cribbed from Zev Kronenberg       #
######################################

CC=gcc
GIT_VERSION := $(shell git describe --abbrev=4 --dirty --always)
FLAGS= -Wall -fopenmp -DVERSION=\"$(GIT_VERSION)\" -std=gnu11 #-fpermissive # -fpermissive is crazy
LD= -lm -lz -lpthread 
INCLUDE= -Isrc -Ihtslib -Ihtslib/htslib -I. -Ilib -Iinclude
LIB=-Lhtslib -lhts

OPT_FLAGS = -O3 -DNDEBUG -finline-functions
DB_FLAGS = -fno-inline -DNOPARALLEL -DTEST
GP_FLAGS = -fno-inline -DNOPARALLEL -DTEST -pg
UR_FLAGS = $(OPT_FLAGS) -DUNROLL

IGAMC_INC= include/igamc_cephes.c

.PHONY: all clean

all: bin lh3sort libhts.a hash_dmp.o fqmarksplit crms bmfsort hash_dmp famstats bam_pr copy

lh3sort:
	cd include/sort && make && cd ../..

bin:
	mkdir -p bin
libhts.a:
	cd htslib && make && cp libhts.a ../
fqmarksplit:
	$(CC) $(FLAGS) $(INCLUDE) $(LIB) $(LD) $(OPT_FLAGS) $(IGAMC_INC) src/fqmarksplit.c  src/khash_dmp_core.c  -o fqmarksplit
	$(CC) $(FLAGS) $(INCLUDE) $(LIB) $(LD) $(DB_FLAGS) $(IGAMC_INC) src/fqmarksplit.c  src/khash_dmp_core.c  -o fqmarksplit_db
	$(CC) $(FLAGS) $(INCLUDE) $(LIB) $(LD) $(DB_FLAGS) $(IGAMC_INC) -DNOPARALLEL  src/khash_dmp_core.c  src/fqmarksplit.c -o fqmarksplit_np
	$(CC) $(FLAGS) $(INCLUDE) $(LIB) $(LD) $(GP_FLAGS) $(IGAMC_INC) src/fqmarksplit.c  src/khash_dmp_core.c  -o fqmarksplit_p
bam_pr:
	$(CC) $(FLAGS) $(INCLUDE) $(LIB) $(LD) $(OPT_FLAGS) src/sam_opts.c src/bam_pr.c -o bam_pr
	$(CC) $(FLAGS) $(INCLUDE) $(LIB) $(LD) $(DB_FLAGS) src/sam_opts.c src/bam_pr.c -o bam_pr_db
	$(CC) $(FLAGS) $(INCLUDE) $(LIB) $(LD) $(PG_FLAGS) src/sam_opts.c src/bam_pr.c -o bam_pr_p
crms:
	$(CC) $(FLAGS) $(INCLUDE) $(LIB) $(LD) $(OPT_FLAGS) src/crms.c include/igamc_cephes.c src/khash_dmp_core.c -o crms
	$(CC) $(FLAGS) $(INCLUDE) $(LIB) $(LD) $(UR_FLAGS) src/crms.c include/igamc_cephes.c src/khash_dmp_core.c -o crms_unroll
	$(CC) $(FLAGS) $(INCLUDE) $(LIB) $(LD) $(DB_FLAGS) src/crms.c include/igamc_cephes.c src/khash_dmp_core.c -o crms_db
	$(CC) $(FLAGS) $(INCLUDE) $(LIB) $(LD) $(GP_FLAGS) src/crms.c include/igamc_cephes.c src/khash_dmp_core.c -o crms_p
bmfsort:
	$(CC) $(FLAGS) $(INCLUDE) $(LIB) $(LD) $(DB_FLAGS) src/bmfsort.c -o bmfsort_db
	$(CC) $(FLAGS) $(INCLUDE) $(LIB) $(LD) $(GP_FLAGS) src/bmfsort.c -o bmfsort_p
	$(CC) $(FLAGS) $(INCLUDE) $(LIB) $(LD) $(OPT_FLAGS) src/bmfsort.c -o bmfsort
hash_dmp.o:
	$(CC) $(FLAGS) $(INCLUDE) $(LIB) $(LD) $(OPT_FLAGS) -fPIC src/khash_dmp_core.c -c -o hash_dmp.o
hash_dmp:
	$(CC) $(FLAGS) $(INCLUDE) $(LIB) $(LD) $(DB_FLAGS) src/khash_dmp_main.c  src/khash_dmp_core.c  include/igamc_cephes.c -o hash_dmp_db
	$(CC) $(FLAGS) $(INCLUDE) $(LIB) $(LD) $(OPT_FLAGS) src/khash_dmp_main.c  src/khash_dmp_core.c  include/igamc_cephes.c -o hash_dmp
	$(CC) $(FLAGS) $(INCLUDE) $(LIB) $(LD) $(UR_FLAGS) src/khash_dmp_main.c  src/khash_dmp_core.c  include/igamc_cephes.c -o hash_dmp_unroll
	$(CC) $(FLAGS) $(INCLUDE) $(LIB) $(LD) $(GP_FLAGS) src/khash_dmp_main.c  src/khash_dmp_core.c  include/igamc_cephes.c -o hash_dmp_p
	#$(CC) $(FLAGS) $(INCLUDE) $(LIB) $(LD) $(DB_FLAGS) src/khash_dmp_main.c  src/khash_dmp_core.c  include/igamc_cephes.c -o hash_dmp_db
	#$(CC) $(FLAGS) $(INCLUDE) $(LIB) $(LD) $(OPT_FLAGS) src/khash_dmp_main.c  src/khash_dmp_core.c  include/igamc_cephes.c -o hash_dmp
	#$(CC) $(FLAGS) $(INCLUDE) $(LIB) $(LD) $(GP_FLAGS) src/khash_dmp_main.c  src/khash_dmp_core.c  include/igamc_cephes.c -o hash_dmp_p
dmp:
	$(CC) $(FLAGS) $(INCLUDE) $(LIB) $(LD) $(OPT_FLAGS) src/dmp.c include/igamc_cephes.c -o dmp
famstats:
	$(CC) $(FLAGS) $(INCLUDE) $(LIB) $(LD) $(OPT_FLAGS) src/famstats.c -o famstats
copy:
	mv hash_dmp crms crms_db crms_p bmfsort bmfsort_db bmfsort_p include/sort/lh3sort fqmarksplit fqmarksplit_db fqmarksplit_p fqmarksplit_np famstats hash_dmp_p hash_dmp_db bam_pr bam_pr_db bam_pr_p bin/


clean:
	rm *.a && rm bin/* && cd htslib && make clean


