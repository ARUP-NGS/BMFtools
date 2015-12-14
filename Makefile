######################################
# Makefile written by Daniel Baker   #
#     d.nephi.baker@gmail.co m       #
#  Cribbed from Zev Kronenberg       #
######################################

CC=gcc
GIT_VERSION := $(shell git describe --abbrev=4 --dirty --always)
FLAGS= -Wall -fopenmp -DVERSION=\"$(GIT_VERSION)\" -std=gnu11
LD= -lm -lz
INCLUDE= -Isrc -Ihtslib -Ihtslib/htslib -I. -Ilib -Iinclude -Isamtools
LIB=

OPT_FLAGS = -O3 -DNDEBUG -flto -fivopts -Wno-unused-function -Wno-unused-variable -Wno-strict-aliasing 
#O2_FLAGS = -O2 -DNDEBUG -finline-functions
DB_FLAGS = -fno-inline -Wno-unused-function -Wno-strict-aliasing
PG_FLAGS = -Wno-unused-function -pg -DNDEBUG -O2 -fno-inline -DPUTC -Wno-strict-aliasing
UR_FLAGS = $(OPT_FLAGS) -DUNROLL

IGAMC_INC = include/igamc_cephes.c

BMF_SRC = htslib/sam.c src/bmf_main.c include/sam_opts.c src/crms.c include/igamc_cephes.c src/hash_dmp_core.c \
		  src/fqmarksplit.c src/bam_rsq.c src/famstats.c src/bmf_vetter.c lib/bed_util.c include/bedidx.c \
		  libhts.a src/bmfsort.c src/err_calc.c lib/io_util.c lib/pair_util.c lib/nix_resource.c \
		  lib/kingfisher.c lib/bam_util.c
		  

.PHONY: all clean

all: bin libhts.a bmftools copy

bin:
	mkdir -p bin
libhts.a:
	cd htslib && make && cp libhts.a ../
bmftools:
	$(CC) $(FLAGS) $(INCLUDE) $(LIB) $(LD) $(OPT_FLAGS) $(BMF_SRC) -o bmftools
	$(CC) $(FLAGS) $(INCLUDE) $(LIB) $(LD) $(DB_FLAGS) $(BMF_SRC) -o bmftools_db
	$(CC) $(FLAGS) $(INCLUDE) $(LIB) $(LD) $(PG_FLAGS) $(BMF_SRC) -o bmftools_p
copy:
	mv bmftools bmftools_p bmftools_db bin/
	#mv bmftools  bin/


clean:
	rm -f *.a && rm -f *.o && rm bin/* # && cd htslib && make clean && cd ../samtools && make clean


