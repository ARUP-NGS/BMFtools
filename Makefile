######################################
# Makefile written by Daniel Baker   #
#     d.nephi.baker@gmail.co m       #
#  Cribbed from Zev Kronenberg       #
######################################

CC=gcc
GIT_VERSION := $(shell git describe --abbrev=4 --dirty --always)
FLAGS= -Wall -fopenmp -DVERSION=\"$(GIT_VERSION)\" -std=gnu99
LD= -lm -lz
INCLUDE= -Ihtslib -I.
LIB=

OPT_FLAGS = -finline-functions -O3 -DNDEBUG -flto -fivopts -Wno-unused-function -Wno-unused-variable -Wno-strict-aliasing 
DB_FLAGS = -Wno-unused-function -Wno-strict-aliasing -fno-inline
PG_FLAGS = -fno-inline -Wno-unused-function -pg -DNDEBUG -O2 -Wno-strict-aliasing
UR_FLAGS = $(OPT_FLAGS) -DUNROLL

IGAMC_INC = include/igamc_cephes.c

OBJS = htslib/sam.o include/sam_opts.o src/bmf_dmp.o include/igamc_cephes.o src/bmf_hashdmp.o \
		  src/bmf_sdmp.o src/bmf_rsq.o src/bmf_famstats.o src/bmf_vetter.o dlib/bed_util.o include/bedidx.o \
		  src/bmf_sort.o src/bmf_err.o dlib/io_util.o dlib/nix_util.o \
		  lib/kingfisher.o dlib/bam_util.o src/bmf_mark_unclipped.o src/bmf_cap.o lib/mseq.o lib/splitter.o

BMF_SRC = htslib/sam.c src/bmf_main.c include/sam_opts.c src/bmf_dmp.c include/igamc_cephes.c src/bmf_hashdmp.c \
		  src/bmf_sdmp.c src/bmf_rsq.c src/bmf_famstats.c src/bmf_vetter.c dlib/bed_util.c include/bedidx.c \
		  libhts.a src/bmf_sort.c src/bmf_err.c dlib/io_util.c dlib/nix_util.c \
		  lib/kingfisher.c dlib/bam_util.c src/bmf_mark_unclipped.c src/bmf_cap.c lib/mseq.c lib/splitter.c
		  

.PHONY: all clean

all: bin libhts.a bmftools copy

%.o: %.c
	$(CC) -c $(FLAGS) $(INCLUDE) $(LIB) $(LD) $(OPT_FLAGS) $< -o $@

bin:
	mkdir -p bin
libhts.a:
	cd htslib && make && cp libhts.a ../
bmftools_db:
	$(CC) $(FLAGS) $(INCLUDE) $(LIB) $(LD) $(DB_FLAGS) $(BMF_SRC) -o bmftools_db
bmftools_p: bmftools_db
	$(CC) $(FLAGS) $(INCLUDE) $(LIB) $(LD) $(PG_FLAGS) $(BMF_SRC) -o bmftools_p
bmftools: $(OBJS)
	$(CC) $(FLAGS) $(INCLUDE) $(LIB) $(LD) $(OPT_FLAGS) $(OBJS) libhts.a src/bmf_main.c -o bmftools
copy: bmftools
	mv bmftools bin/ #mv bmftools bmftools_p bmftools_db bin/


clean:
	rm -f *.a && rm -f *.o && rm -f bin/* && rm -f src/*.o && rm -f dlib/*.o && \
		rm -f include/*.o && rm -f lib/*.o && cd htslib && make clean


