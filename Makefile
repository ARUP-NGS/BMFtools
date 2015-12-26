######################################
# Makefile written by Daniel Baker   #
#     d.nephi.baker@gmail.com       #
######################################

CC=gcc
GIT_VERSION := $(shell git describe --abbrev=4 --dirty --always)
FLAGS= -Wall -fopenmp -DVERSION=\"$(GIT_VERSION)\" -std=gnu99
LD= -lm -lz
INCLUDE= -Ihtslib -Iinclude -I.
LIB=
INSTALL=/usr/bin/install -c

prefix = /usr/local
bindir = $(prefix)/bin
binprefix =

OPT_FLAGS = -finline-functions -O3 -DNDEBUG -flto -fivopts -Wno-unused-function -Wno-unused-variable -Wno-strict-aliasing 
DB_FLAGS = -Wno-unused-function -Wno-strict-aliasing -fno-inline
PG_FLAGS = -fno-inline -Wno-unused-function -pg -DNDEBUG -O2 -Wno-strict-aliasing
UR_FLAGS = $(OPT_FLAGS) -DUNROLL

IGAMC_INC = include/igamc_cephes.c

OBJS = htslib/sam.o include/sam_opts.o src/bmf_dmp.o include/igamc_cephes.o src/bmf_hashdmp.o \
		  src/bmf_sdmp.o src/bmf_rsq.o src/bmf_famstats.o src/bmf_vetter.o dlib/bed_util.o include/bedidx.o \
		  src/bmf_sort.o src/bmf_err.o dlib/io_util.o dlib/nix_util.o \
		  lib/kingfisher.o dlib/bam_util.o src/bmf_mark_unclipped.o src/bmf_cap.o lib/mseq.o lib/splitter.o \
		  src/bmf_main.o

# In case you want to make a debug or profile build without changing the .o/.c rules.
BMF_SRC = $(OBJS:.o=.c) libhts.a

.PHONY: all clean install

all: libhts.a bmftools
install: all
	$(INSTALL) bmftools $(bindir)/$(binprefix)bmftools
	$(INSTALL) bmftools_p $(bindir)/$(binprefix)bmftools_p
	$(INSTALL) bmftools_db $(bindir)/$(binprefix)bmftools_db

%.o: %.c
	$(CC) -c $(FLAGS) $(INCLUDE) $(LIB) $(LD) $(OPT_FLAGS) $< -o $@

libhts.a:
	cd htslib && make && cp libhts.a ../
bmftools_db: libhts.a
	$(CC) $(FLAGS) $(INCLUDE) $(LIB) $(LD) $(DB_FLAGS) $(BMF_SRC) -o bmftools_db
bmftools_p: bmftools_db
	$(CC) $(FLAGS) $(INCLUDE) $(LIB) $(LD) $(PG_FLAGS) $(BMF_SRC) -o bmftools_p
bmftools: $(OBJS) bmftools_p
	$(CC) $(FLAGS) $(INCLUDE) $(LIB) $(LD) $(OPT_FLAGS) $(OBJS) libhts.a -o bmftools


clean:
	rm -f *.a && rm -f *.o && rm -f bmftools* && rm -f src/*.o && rm -f dlib/*.o && \
		rm -f include/*.o && rm -f lib/*.o && cd htslib && make clean


