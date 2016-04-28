######################################
# Makefile written by Daniel Baker   #
#     d.nephi.baker@gmail.com       #
######################################


CXXSTD=c++11
CSTD=gnu99
CC=g++
GIT_VERSION := $(shell git describe --abbrev=4 --dirty --always)
CFLAGS= -Wunreachable-code -Wall -fopenmp -DBMF_VERSION=\"$(GIT_VERSION)\" -std=$(CSTD) -fno-builtin-gamma -pedantic
FLAGS= -Wunreachable-code -Wall -fopenmp -DBMF_VERSION=\"$(GIT_VERSION)\" -std=$(CXXSTD) -fno-builtin-gamma -pedantic
LD= -lm -lz -lpthread
INCLUDE= -Ihtslib -Iinclude -I.
LIB=
INSTALL=/usr/bin/install -c
THREADS=12
GENOME_PATH=/mounts/genome/human_g1k_v37.fasta

prefix = /usr/local
bindir = $(prefix)/bin
binprefix =

OPT_FLAGS = -finline-functions -O3 -DNDEBUG -flto -fivopts -Wno-unused-function -Wno-strict-aliasing -fno-builtin-gamma
DB_FLAGS = -Wno-unused-function -Wno-strict-aliasing -pedantic -fno-builtin-gamma -fno-inline
PG_FLAGS = -Wno-unused-function -pg -DNDEBUG -O3 -Wno-strict-aliasing -fno-builtin-gamma -fno-inline

DLIB_SRC = dlib/cstr_util.c dlib/math_util.c dlib/vcf_util.c dlib/io_util.c dlib/bam_util.c dlib/nix_util.c \
		   dlib/bed_util.c dlib/misc_util.c

SOURCES = include/sam_opts.c src/bmf_dmp.c include/igamc_cephes.c src/bmf_hashdmp.c \
          src/bmf_sdmp.c src/bmf_rsq.c src/bmf_famstats.c include/bedidx.c \
          src/bmf_err.c src/bmf_infer.c\
          lib/kingfisher.c src/bmf_mark.c src/bmf_cap.c lib/mseq.c lib/splitter.c \
          src/bmf_main.c src/bmf_target.c src/bmf_depth.c src/bmf_vetter.c src/bmf_sort.c src/bmf_stack.c \
          lib/stack.c src/bmf_filter.c $(DLIB_SRC)

TEST_SOURCES = test/target_test.c test/ucs/ucs_test.c test/tag/array_tag_test.c

TEST_OBJS = $(TEST_SOURCES:.c=.dbo)

P_OBJS = $(SOURCES:.c=.po)
D_OBJS = $(SOURCES:.c=.dbo)
OBJS = $(SOURCES:.c=.o)
DLIB_OBJS = $(DLIB_SRC:.c=.o)


ALL_TESTS=test/ucs/ucs_test marksplit_test hashdmp_test target_test err_test rsq_test
BINS=bmftools bmftools_db bmftools_p

.PHONY: all clean install tests python mostlyclean hashdmp_test err_test update_dlib

all: update_dlib libhts.a tests $(BINS)

install: all
	$(INSTALL) bmftools $(bindir)/$(binprefix)bmftools
	$(INSTALL) bmftools_db $(bindir)/$(binprefix)bmftools_db
	$(INSTALL) bmftools_p $(bindir)/$(binprefix)bmftools_p

%.o: %.cpp
	$(CC) -c $(FLAGS) $(INCLUDE) $(LIB) $(LD) $(OPT_FLAGS) $< -o $@

src/%.o: src/%.cpp cstr_util.o
	$(CC) -c $(FLAGS) $(INCLUDE) $(LIB) $(LD) $(OPT_FLAGS) $(DLIB_OBJS) $< -o $@

%.o: %.c
	gcc -c $(CFLAGS) $(INCLUDE) $(LIB) $(LD) $(OPT_FLAGS) $< -o $@

%.po: %.cpp
	$(CC) -c $(FLAGS) $(INCLUDE) $(LIB) $(LD) $(PG_FLAGS) $< -o $@

%.po: %.c
	gcc -c $(CFLAGS) $(INCLUDE) $(LIB) $(LD) $(PG_FLAGS) $< -o $@

%.dbo: %.cpp
	$(CC) -c $(FLAGS) $(INCLUDE) $(LIB) $(LD) $(DB_FLAGS) $< -o $@

%.dbo: %.c
	gcc -c $(CFLAGS) $(INCLUDE) $(LIB) $(LD) $(DB_FLAGS) $< -o $@


libhts.a:
	+cd htslib && echo "/* Empty config.h */" >> config.h && make -j $(THREADS) && cp libhts.a ../
bmftools_db: $(D_OBJS) libhts.a
	$(CC) $(FLAGS) $(INCLUDE) $(LIB) $(LD) $(DB_FLAGS) $(D_OBJS) libhts.a -o bmftools_db
bmftools_p: $(P_OBJS) libhts.a
	$(CC) $(FLAGS) $(INCLUDE) $(LIB) $(LD) $(PG_FLAGS) $(P_OBJS) libhts.a -o bmftools_p
bmftools: $(OBJS) libhts.a
	$(CC) $(FLAGS) $(INCLUDE) $(LIB) $(LD) $(OPT_FLAGS) $(OBJS) libhts.a -o bmftools
test/ucs/ucs_test: libhts.a $(TEST_OBJS)
	$(CC) $(FLAGS) $(INCLUDE) $(LIB) $(LD) $(DB_FLAGS) test/ucs/ucs_test.dbo libhts.a -o test/ucs/ucs_test
	cd test/ucs && ./ucs_test && cd ./..
tag_test: $(OBJS) $(TEST_OBJS) libhts.a
	$(CC) $(FLAGS) $(DB_FLAGS) $(INCLUDE) $(LIB) $(LD) test/tag/array_tag_test.dbo libhts.a -o ./tag_test && ./tag_test
target_test: $(D_OBJS) $(TEST_OBJS) libhts.a
	$(CC) $(FLAGS) $(DB_FLAGS) $(INCLUDE) $(LIB) $(LD) dlib/bed_util.dbo src/bmf_target.dbo test/target_test.dbo libhts.a -o ./target_test && ./target_test
hashdmp_test: $(BINS)
	cd test/dmp && python hashdmp_test.py && cd ../..
marksplit_test: $(BINS)
	cd test/marksplit && python marksplit_test.py && cd ../..
err_test: $(BINS)
	cd test/err && python err_test.py $(GENOME_PATH) && cd ../..
rsq_test: $(BINS)
	cd test/rsq && python rsq_test.py  && cd ../..

tests: $(BINS) $(ALL_TESTS) test/tag/array_tag_test.dbo
	@echo "Passed all tests!"

python:
	cd CyBMFtools && python setup.py build_ext && python setup.py install && cd ..


clean: mostlyclean
		cd htslib && make clean && cd ..

update_dlib:
	cd dlib && git checkout master && git pull origin master && cd ..

mostlyclean:
	rm -f *.*o && rm -f bmftools* && rm -f src/*.*o && rm -f dlib/*.*o && \
		rm -f include/*.*o && rm -f lib/*.*o && rm -f $(find ./test -name '*o')
