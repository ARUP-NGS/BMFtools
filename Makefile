######################################
# Makefile written by Daniel Baker   #
#     d.nephi.baker@gmail.co m       #
#  Cribbed from Zev Kronenberg       #
######################################

CC=gcc
GIT_VERSION := $(shell git describe --abbrev=4 --dirty --always)
FLAGS= -Wall -fopenmp -DVERSION=\"$(GIT_VERSION)\" -std=gnu99
LD= -lm -lz
INCLUDE= -Isrc -Ihtslib -Ihtslib/htslib -I. -Ilib -Iinclude -Isamtools -Idlib
LIB=

OPT_FLAGS = -finline-functions -O3 -DNDEBUG -flto -fivopts -Wno-unused-function -Wno-unused-variable -Wno-strict-aliasing 
#O2_FLAGS = -O2 -DNDEBUG -finline-functions
DB_FLAGS = -fno-inline -Wno-unused-function -Wno-strict-aliasing
PG_FLAGS = -fno-inline -Wno-unused-function -pg -DNDEBUG -O2 -Wno-strict-aliasing
UR_FLAGS = $(OPT_FLAGS) -DUNROLL

IGAMC_INC = include/igamc_cephes.c

BMF_SRC = htslib/sam.c src/bmf_main.c include/sam_opts.c src/bmf_dmp.c include/igamc_cephes.c src/bmf_hashdmp.c \
		  src/bmf_sdmp.c src/bmf_rsq.c src/bmf_famstats.c src/bmf_vetter.c dlib/bed_util.c include/bedidx.c \
		  libhts.a src/bmf_sort.c src/bmf_err.c dlib/io_util.c dlib/nix_util.c \
		  lib/kingfisher.c dlib/bam_util.c src/bmf_mark_unclipped.c src/bmf_cap.c lib/mseq.c lib/splitter.c
		  

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


