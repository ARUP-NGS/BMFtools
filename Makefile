######################################
# Makefile written by Daniel Baker   #
#     d.nephi.baker@gmail.com       #
######################################

CC=gcc
GIT_VERSION := $(shell git describe --abbrev=4 --dirty --always)
FLAGS= -Wall -fopenmp -DVERSION=\"$(GIT_VERSION)\" -std=gnu99 -pedantic
LD= -lm -lz -lpthread
INCLUDE= -Ihtslib -Iinclude -I.
LIB=
INSTALL=/usr/bin/install -c
THREADS=12

prefix = /usr/local
bindir = $(prefix)/bin
binprefix =

OPT_FLAGS = -finline-functions -O3 -DNDEBUG -flto -fivopts -Wno-unused-function -Wno-unused-variable -Wno-strict-aliasing
DB_FLAGS = -Wno-unused-function -Wno-strict-aliasing
PG_FLAGS = -Wno-unused-function -pg -DNDEBUG -O3 -Wno-strict-aliasing

SOURCES = htslib/sam.c include/sam_opts.c src/bmf_dmp.c include/igamc_cephes.c src/bmf_hashdmp.c \
		  src/bmf_sdmp.c src/bmf_rsq.c src/bmf_famstats.c src/bmf_vetter.c dlib/bed_util.c include/bedidx.c \
		  src/bmf_sort.c src/bmf_err.c dlib/io_util.c dlib/nix_util.c \
		  lib/kingfisher.c dlib/bam_util.c src/bmf_mark_unclipped.c src/bmf_cap.c lib/mseq.c lib/splitter.c \
		  src/bmf_main.c src/bmf_target.c src/bmf_depth.c

TEST_SOURCES = test/target_test.c

TEST_OBJS = $(TEST_SOURCES:.c=.o)

P_OBJS = $(SOURCES:.c=.po)
D_OBJS = $(SOURCES:.c=.dbo)
OBJS = $(SOURCES:.c=.o)

.PHONY: all clean install tests python mostlyclean

all: libhts.a tests bmftools bmftools_db bmftools_p

install: all
	$(INSTALL) bmftools $(bindir)/$(binprefix)bmftools
	$(INSTALL) bmftools_db $(bindir)/$(binprefix)bmftools_db
	$(INSTALL) bmftools_p $(bindir)/$(binprefix)bmftools_p

%.o: %.c
	$(CC) -c $(FLAGS) $(INCLUDE) $(LIB) $(LD) $(OPT_FLAGS) $< -o $@

%.po: %.c
	$(CC) -c $(FLAGS) $(INCLUDE) $(LIB) $(LD) $(PG_FLAGS) $< -o $@

%.dbo: %.c
	$(CC) -c $(FLAGS) $(INCLUDE) $(LIB) $(LD) $(DB_FLAGS) $< -o $@

libhts.a:
	cd htslib && make -j $(THREADS) && cp libhts.a ../
bmftools_db: $(D_OBJS) libhts.a
	$(CC) $(FLAGS) $(INCLUDE) $(LIB) $(LD) $(DB_FLAGS) $(D_OBJS) libhts.a -o bmftools_db
bmftools_p: $(P_OBJS) libhts.a
	$(CC) $(FLAGS) $(INCLUDE) $(LIB) $(LD) $(PG_FLAGS) $(P_OBJS) libhts.a -o bmftools_p
bmftools: $(OBJS) libhts.a
	$(CC) $(FLAGS) $(INCLUDE) $(LIB) $(LD) $(OPT_FLAGS) $(OBJS) libhts.a -o bmftools
test/ucs/ucs_test: libhts.a
	$(CC) $(FLAGS) $(INCLUDE) $(LIB) $(LD) $(DB_FLAGS) test/ucs/ucs_test.c libhts.a -o test/ucs/ucs_test


tests: $(TEST_OBJS) dlib/bed_util.o src/bmf_target.o bmftools bmftools_db bmftools_p test/ucs/ucs_test
	$(CC) $(FLAGS) $(DB_FLAGS) $(INCLUDE) $(LIB) $(LD) dlib/bed_util.o src/bmf_target.o test/target_test.c libhts.a -o ./target_test && ./target_test
	cd test/ucs && ./ucs_test && cd ../..
	cd test/dmp && python hashdmp_test.py && cd ../..
	cd test/marksplit && python marksplit_test.py && cd ../..

python:
	cd CyBMFtools && python setup.py build_ext && python setup.py install && cd ..



clean:
	rm -f *.a && rm -f *.o && rm -f bmftools* && rm -f src/*.o && rm -f dlib/*.o && \
		rm -f include/*.o && rm -f lib/*.o && cd htslib && make clean

mostlyclean:
	rm -f *.*o && rm -f bmftools* && rm -f src/*.*o && rm -f dlib/*.*o && \
		rm -f include/*.*o && rm -f lib/*.*o && rm -f test/*o
