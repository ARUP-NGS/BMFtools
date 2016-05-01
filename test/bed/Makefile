CXX=g++
CXXSTD=c++11
CC=g++
SRC=bed_test.cpp ../../include/bedidx.c ../../dlib/bed_util.cpp
OBJtmp=$(SRC:.cpp=.o)
OBJ=$(OBJtmp:.c=.co)
INCLUDE=-I ../../
FLAGS= -Wunreachable-code -Wall -fopenmp -DBMF_VERSION=\"$(GIT_VERSION)\" -std=$(CXXSTD) -fno-builtin-gamma -pedantic
LD= -lm -lz -lpthread
INCLUDE= -I../.. -I../..htslib -I../../include -I.
LIB=

.PHONY: all clean
all: libhts.a bed_test

%.o: %.cpp
	$(CXX) -c $(FLAGS) $(INCLUDE) $(LIB) $(LD) $(OPT_FLAGS) $< -o $@
%.co: %.c
	$(CC) -c $(FLAGS) $(INCLUDE) $(LIB) $(LD) $(OPT_FLAGS) $< -o $@

libhts.a:
	cd ../../htslib && make -j14 && cd ../test/bed && ln -s ../../htslib/libhts.a libhts.a

bed_test: $(OBJ) ../../libhts.a
	$(CXX) $(FLAGS) $(INCLUDE) $(LIB) $(LD) $(OBJ) ../../libhts.a -o bed_test

clean:
	rm -f ./*.o libhts.a && cd ../../htslib/ && make clean && cd ../test/bed

